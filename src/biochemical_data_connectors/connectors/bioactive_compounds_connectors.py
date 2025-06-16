import time
import requests
import logging
import concurrent.futures
from abc import ABC, abstractmethod
from functools import partial
from typing import List, Dict, Optional, Any

import pubchempy as pcp
from chembl_webresource_client.new_client import new_client

from biochemical_data_connectors.constants import RestApiEndpoints
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.utils.iter_utils import batch_iterable
from biochemical_data_connectors.utils.api.mappings import uniprot_to_gene_id_mapping
from biochemical_data_connectors.utils.api.pubchem_api import PubChemAPIClient, get_compounds_in_batches


class BaseBioactivesConnector(ABC):
    """
    Abstract base class for extracting bioactive compounds from a data source.

    Attributes
    ----------
    _bioactivity_measures : List[str]
        A prioritized list of bioactivity measurement types to filter on
        (e.g., ['Kd', 'Ki', 'IC50']).
    _bioactivity_threshold : float, optional
        The maximum potency value (in nM) to consider a compound bioactive.
    _logger : logging.Logger
        A logger instance for logging messages.

    Methods
    -------
    get_bioactive_compounds(target: str) -> List[str]
        Abstract method to return a list of `BioactiveCompound` objects for a target UniProt ID
        identifier.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,
        logger: Optional[logging.Logger] = None
    ):
        self._bioactivity_measures = bioactivity_measures
        self._bioactivity_threshold = bioactivity_threshold
        self._logger = logger if logger else logging.getLogger(__name__)

    @abstractmethod
    def get_bioactive_compounds(self, target_uniprot_id: str) -> List[BioactiveCompound]:
        """
        Retrieve a list of canonical SMILES for bioactive compounds for a given target.

        Parameters
        ----------
        target_uniprot_id : str
            The target identifier (UniProt accession, e.g. "P00533").

        Returns
        -------
        List[BioactiveCompound]
            A list of structured BioactiveCompound objects.
        """
        pass


class ChEMBLBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts bioactive compounds for a given target from ChEMBL using a UniProt accession.

    Attributes
    ----------
    _client : object, optional
        A ChEMBL client instance used for making API requests. If None, the default client
        from `chembl_webresource_client.new_client` will be used.

    Methods
    -------
    get_bioactive_compounds(target: str) -> List[str]
        Returns a list of `BioactiveCompound` objects for a target UniProt ID
        identifier.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        client = None,
        bioactivity_threshold: Optional[float] = None, # In nM (e.g. 1000 nM threshold to filter for compounds with Kd <= 1 µM)
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(bioactivity_measures, bioactivity_threshold, logger)
        self._client = client if client else new_client

    def get_bioactive_compounds(self, target_uniprot_id: str) -> List[BioactiveCompound]:
        """
        Retrieve bioactive compound objects for a given target from ChEMBL.

        This method queries the ChEMBL activity API, fetching full records
        for compounds that match the target and bioactivity criteria, and
        returns them as a list of structured BioactiveCompound objects.

        Parameters
        ----------
        target_uniprot_id : str
            The UniProt accession for the target.

        Returns
        -------
        List[BioactiveCompound]
            A list of BioactiveCompound objects meeting the criteria.
        """
        # 1) Search for the target by UniProt ID and retrieve the first matching result
        target_results = self._client.target.filter(target_components__accession=target_uniprot_id)
        if not target_results:
            self._logger.error(f"No matching target found for UniProt ID {target_uniprot_id}")
            return []

        target_id = target_results[0]['target_chembl_id']

        # 2) Build the base parameters for ChEMBL REST API.
        #    Filter activities based on the specified standard type (e.g. IC50)
        params = {
            "target_chembl_id": target_id,
            "standard_type": self._bioactivity_measures,
            "standard_units": "nM"
        }

        # 2.1) Add threshold value to params if set
        if self._bioactivity_threshold is not None:
            # standard_value__lte is the ChEMBL REST parameter for “≤”
            params["standard_value__lte"] = self._bioactivity_threshold

        # 2.2) Paginate REST API requests
        limit: int = 1000
        offset: int = 0
        chembl_activity_url: str = RestApiEndpoints.CHEMBL_ACTIVITY.url()
        bioactive_compounds: List[BioactiveCompound] = []
        chembl_start = time.time()

        while True:
            page_params = {
                **params,
                "limit": limit,
                "offset": offset
            }

            chembl_activity_request = requests.get(
                chembl_activity_url,
                params=page_params,
                timeout=15
            )

            chembl_activity_request.raise_for_status()
            activity_data = chembl_activity_request.json()

            records = activity_data.get('activities', [])
            if not records:
                break

            for record in records:
                structures = record.get('molecule_structures')
                properties = record.get('molecule_properties')

                if not all([
                    record.get('molecule_chembl_id'),
                    record.get('canonical_smiles'),
                    record.get('standard_value')
                ]):
                    continue

                compound_obj = BioactiveCompound(
                    source_db="ChEMBL",
                    source_id=record['molecule_chembl_id'],
                    smiles=record['canonical_smiles'],
                    activity_type=record['standard_type'],
                    activity_value=float(record['standard_value']),
                    source_inchikey=structures['standard_inchi_key'] if structures['standard_inchi_key'] else None,
                    iupac_name=structures.get('iupac_name') if structures.get('iupac_name') else None,
                    molecular_formula=properties.get('full_molformula') if properties and properties.get(
                        'full_molformula') else None,
                    molecular_weight=float(properties.get('mw_freebase')) if properties and properties.get(
                        'mw_freebase') else None,
                    raw_data=record
                )
                bioactive_compounds.append(compound_obj)

            offset += limit

        chembl_end = time.time()
        self._logger.info(f'ChEMBL total query time: {round(chembl_end - chembl_start)} seconds')

        return bioactive_compounds


class PubChemBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts bioactive compounds for a given target from PubChem using a UniProt accession.

    This connector orchestrates a multi-step process to query PubChem,
    retrieve all relevant compounds and their bioactivity data, and formats
    them into standardized `BioactiveCompound` objects.

    Methods
    -------
    get_bioactive_compounds(target: str) -> List[str]
        Retrieves canonical SMILES for compounds from PubChem for the given UniProt target.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None, # In nM (e.g. 1000 nM threshold to filter for compounds with Kd <= 1 µM)
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(bioactivity_measures, bioactivity_threshold, logger)
        self._api_client = PubChemAPIClient(logger=self._logger)

    def get_bioactive_compounds(self, target_uniprot_id: str) -> List[BioactiveCompound]:
        """
        Retrieve canonical SMILES for compounds for a given target from PubChem.
        The target is provided as a UniProt accession (e.g. "P00533").

        This method performs the following steps:

        1. Maps the UniProt accession to an NCBI GeneID.
        2. Uses the GeneID to query PubChem’s BioAssay API for assay IDs (AIDs).
        3. For each assay, extracts the active compound IDs (CIDs).
        4. Retrieves full `pubchempy.Compound` compound details for all CIDs.
        5. Creates a standardized BioactiveCompound object for each compound
           that has a valid potency score.
        6. Optionally filters the final list based on the potency threshold.

        Parameters
        ----------
        target_uniprot_id : str
            The UniProt accession for the target (e.g., "P00533").

        Returns
        -------
        List[BioactiveCompound]
            A list of standardized BioactiveCompound objects
        """
        # 1) Map the UniProt accession to an NCBI GeneID.
        target_gene_id = self._lookup_target_gene_id(target_uniprot_id)
        if not target_gene_id:
            self._logger.error(f"Could not determine GeneID for target '{target_uniprot_id}'.")
            return []

        # 2) Query the BioAssay API to get assay IDs (AIDs) for the target GeneID.
        try:
            aid_list = self._api_client.get_active_aids(target_gene_id)
        except Exception as e:
            self._logger.error(f"Error retrieving assay IDs for GeneID {target_gene_id}: {e}")
            return []

        # 3) For each assay, retrieve active compound IDs (CIDs) and aggregate them.
        #    Create thread pool using Python’s `ThreadPoolExecutor` to issue multiple API calls concurrently in batches
        cids_api_start: float = time.time()
        active_cids = set()

        # Create a new partial function with `logger` argument fixed. This allows us to pass a fixed `logger` argument
        # to the `get_active_cids_wrapper()` function when it is mapped to each AID element in `aid_list` via
        # `concurrent.futures.ThreadPoolExecutor.map()`
        get_active_cids_partial = partial(self._api_client.get_active_cids, logger=self._logger)

        with concurrent.futures.ThreadPoolExecutor(max_workers=9) as executor:
            # Map and apply partial function of `cids_for_aid_wrapper()` to every element in `aid_list` concurrently
            results = list(executor.map(get_active_cids_partial, aid_list))

            for cids in results:
                active_cids.update(cids)

        cids_api_end: float = time.time()
        self._logger.info(f'PubChem CID total API query time: {round(cids_api_end - cids_api_start)} seconds')

        if not active_cids:
            self._logger.error(f"No active compounds found for GeneID {target_gene_id}.")
            return []

        # 4) Retrieve full `pubchempy.Compound` objects for all CIDs.
        pubchempy_compound_api_start: float = time.time()
        pubchempy_compounds = get_compounds_in_batches(cids=list(active_cids), logger=self._logger)
        pubchempy_compound_api_end: float = time.time()
        self._logger.info(f'PubChem bioactive compounds from CIDs total API query time: '
                          f'{round(pubchempy_compound_api_end - pubchempy_compound_api_start)} seconds')

        # 5) Fetch potencies for ALL retrieved compounds.
        self._logger.info(f"Fetching potencies for {len(pubchempy_compounds)} compounds...")
        potencies_api_start: float = time.time()
        cid_to_potency_map = {}

        # Create a new partial function with `target_gene_id` and `logger` argument fixed. As before, this allows
        # us to pass these fixed arguments to `self._get_compound_bioassay_data()` when it is mapped to each
        # compound element in the batched `bioactive_compounds` iterable via
        # `concurrent.futures.ThreadPoolExecutor.map()`
        get_compound_bioassay_data_partial = partial(
            self._api_client.get_compound_bioassay_data,
            target_gene_id=target_gene_id,
            bioactivity_measures=self._bioactivity_measures,
            logger=self._logger
        )
        for compound_batch in batch_iterable(iterable=pubchempy_compounds):
            # Process the current `bioactive_compounds` batch concurrently using a thread pool
            with (concurrent.futures.ThreadPoolExecutor(max_workers=9) as executor):
                # Map and apply partial function of `self._get_compound_bioassay_data()` to every element in
                # current `bioactive_compounds` batch concurrently
                batch_bioassay_data = list(
                    executor.map(
                        get_compound_bioassay_data_partial,
                        compound_batch
                    )
                )
                for compound, bioassay_data in zip(compound_batch, batch_bioassay_data):
                    if bioassay_data:
                        cid_to_potency_map[compound.cid] = bioassay_data

        potencies_api_end: float = time.time()
        self._logger.info(f"PubChem bioactive compound potencies total API query time: "
                          f"{round(potencies_api_end - potencies_api_start)} seconds\n"
                          f"Found potency data for {len(cid_to_potency_map)} compounds.")

        # 6) Create complete list of `BioactiveCompound` objects
        all_bioactives: List[BioactiveCompound] = []
        for pubchempy_compound in pubchempy_compounds:
            compound_bioassay_data: Dict = cid_to_potency_map.get(pubchempy_compound.cid)

            if compound_bioassay_data is None:
                self._logger.debug(f"Skipping compound CID {pubchempy_compound.cid} due to missing potency data.")
                continue

            compound_obj = self._create_bioactive_compound(
                pubchempy_compound=pubchempy_compound,
                bioassay_data=compound_bioassay_data
            )
            if compound_obj:
                all_bioactives.append(compound_obj)

        # 7) Filter final list of `BioactiveCompound` objects by potency if threshold is provided.
        if self._bioactivity_threshold is not None:
            self._logger.info(f"Filtering {len(all_bioactives)} compounds with threshold: "
                              f"<= {self._bioactivity_threshold} nM")
            filtered_bioactives: List[BioactiveCompound] = [
                compound for compound in all_bioactives if compound.activity_value <= self._bioactivity_threshold
            ]
            self._logger.info(f"Found {len(filtered_bioactives)} compounds after filtering.")

            return filtered_bioactives
        else:
            return all_bioactives

    @staticmethod
    def _create_bioactive_compound(
        pubchempy_compound: pcp.Compound,
        bioassay_data: Dict[str, Any]
    ) -> Optional[BioactiveCompound]:
        """
        Helper to convert a `pubchempy.Compound` to a `BioactiveCompound`.

        This method safely extracts attributes from the source object and uses
        them to instantiate the standardized `BioactiveCompound` dataclass.

        Parameters
        ----------
        pubchempy_compound : pcp.Compound
            The source object from the `pubchempy` library.
        bioassay_data : Dict[str, Any]
            The dictionary of pre-fetched bioassay data.

        Returns
        -------
        Optional[BioactiveCompound]
            A populated `BioactiveCompound` object, or None if essential
            information like SMILES is missing.
        """
        if not getattr(pubchempy_compound, 'canonical_smiles', None):
            return None

        return BioactiveCompound(
            source_db='PubChem',
            source_id=pubchempy_compound.cid,
            smiles=pubchempy_compound.canonical_smiles,
            activity_type=bioassay_data['activity_type'],
            activity_value=bioassay_data['best_value'],
            source_inchikey=pubchempy_compound.inchikey if pubchempy_compound.inchikey else None,
            iupac_name=getattr(pubchempy_compound, 'iupac_name', None),
            molecular_formula=getattr(pubchempy_compound, 'molecular_formula', None),
            molecular_weight=float(pubchempy_compound.molecular_weight) if getattr(pubchempy_compound,
                                                                                   'molecular_weight', None) else None,
            n_measurements=bioassay_data["n_measurements"],
            mean_activity=bioassay_data["mean_value"],
            median_activity=bioassay_data["median_value"],
            std_dev_activity=bioassay_data["std_dev_value"],
            raw_data=pubchempy_compound
        )

    @staticmethod
    def _lookup_target_gene_id(target: str) -> Optional[str]:
        """
        Look up the target gene identifier (GeneID) for the given UniProt accession by
        using the UniProt ID mapping API.

        Parameters
        ----------
        target : str
            The UniProt accession (e.g., "P00533").

        Returns
        -------
        Optional[str]
            The corresponding NCBI GeneID if found, otherwise None.
        """
        return uniprot_to_gene_id_mapping(target)
