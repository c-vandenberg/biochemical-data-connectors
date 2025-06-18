import os
import time
import json
import pickle
import logging
import statistics
import concurrent.futures
from abc import ABC, abstractmethod
from datetime import datetime
from functools import partial
from collections import defaultdict
from typing import List, Dict, Optional, Any, Callable

import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from biochemical_data_connectors.constants import RestApiEndpoints, CONVERSION_FACTORS_TO_NM
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.utils.iter_utils import batch_iterable
from biochemical_data_connectors.utils.api.chembl_api import ChEMBLAPIClient
from biochemical_data_connectors.utils.api.pubchem_api import PubChemAPIClient, get_compounds_in_batches
from biochemical_data_connectors.utils.api.mappings import uniprot_to_gene_id_mapping

CHEMBL_INVALID_DATA_COMMENT = 'OUTSIDE TYPICAL RANGE'


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
        cache_dir: str = './data/cache',
        logger: Optional[logging.Logger] = None
    ):
        self._bioactivity_measures = bioactivity_measures
        self._bioactivity_threshold = bioactivity_threshold
        self._cache_dir = cache_dir
        self._logger = logger if logger else logging.getLogger(__name__)

    @abstractmethod
    def get_bioactive_compounds(self, target_uniprot_id: str, force_refresh: bool = False) -> List[BioactiveCompound]:
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

    def _get_cached_or_fetch(
        self,
        cache_file_path: str,
        fetch_function: Callable[[], Any],
        data_type: str = 'Bioactive',
        use_pickle: bool = False,
        force_refresh: bool = False
    ) -> Any:
        """
        Generic method to retrieve data from a cache file or by executing a
        data fetching function if the cache does not exist or a refresh is forced.

        Parameters
        ----------
        cache_file_path : str
            The path to the cache file.
        fetch_function : Callable[[], Any]
            A no-argument function that will be called to get fresh data if needed.
        force_refresh : bool, optional
            If True, ignores the cache and always executes the fetch_function.
        use_pickle : bool, optional
            If True, use the binary pickle format for serialization. If False
            (default), use JSON.

        Returns
        -------
        Any
            The data from the cache or the fetch_function.
        """
        cache_is_valid = False
        # 1) If manual `--force-refresh` flag is True, and cache file exists, load data from
        #    cache file.
        if not force_refresh and os.path.exists(cache_file_path):
            if use_pickle:
                with open(cache_file_path, 'rb') as f:  # 'rb' for read binary
                    cache_data = pickle.load(f)
            else:
                with open(cache_file_path, 'r') as f:
                    cache_data = json.load(f)

            self._logger.info(f"Found valid cache at {cache_file_path}. Loading {data_type} data from file.")
            data = cache_data["data"]
            cache_is_valid = True

        # 2) If cache file doesn't exist or `force_refresh == True`, call the provided fetch_function.
        if not cache_is_valid:
            self._logger.info(f"Fetching fresh {data_type} data from API...")

            data = fetch_function()

            # 3) Save the new results to the cache file with a current timestamp.
            if data:
                cache_content = {
                    "timestamp": datetime.now().isoformat(),
                    "data": data
                }
                os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
                if use_pickle:
                    with open(cache_file_path, 'wb') as f:
                        pickle.dump(cache_content, f)
                else:
                    with open(cache_file_path, 'w') as f:
                        json.dump(cache_content, f, indent=4)

                self._logger.info(f"Saved {len(data)} {data_type} items to cache file: {cache_file_path}")

        return data


class ChEMBLBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts bioactive compounds from ChEMBL using a target's UniProt accession.

    Attributes
    ----------
    _chembl_webresource_client : object
        A client for the high-level ChEMBL API, used for target lookups.
    _chembl_api_client : ChEMBLAPIClient
        A client for the low-level ChEMBL REST API, used for activity fetching.

    Methods
    -------
    get_bioactive_compounds(target: str) -> List[BioactiveCompound]
        Returns a list of `BioactiveCompound` objects for a target UniProt ID
        identifier.
    """

    def __init__(
        self,
        bioactivity_measures: List[str],
        core_chembl_client=None,
        bioactivity_threshold: Optional[float] = None,
        # In nM (e.g. 1000 nM threshold to filter for compounds with Kd <= 1 µM)
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(bioactivity_measures, bioactivity_threshold, logger)
        self._chembl_webresource_client = core_chembl_client if core_chembl_client else new_client
        self._chembl_api_client: ChEMBLAPIClient = ChEMBLAPIClient(logger=self._logger)

    def get_bioactive_compounds(self, target_uniprot_id: str, force_refresh: bool = False) -> List[BioactiveCompound]:
        """
        Retrieve bioassay data for bioactive compounds from ChEMBL using a target's UniProt accession.

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
        target_results = self._chembl_webresource_client.target.filter(target_components__accession=target_uniprot_id)
        if not target_results:
            self._logger.error(f"No matching target found for UniProt ID {target_uniprot_id}")
            return []

        target_chembl_id = target_results[0]['target_chembl_id']

        os.makedirs(self._cache_dir, exist_ok=True)
        pubchem_acitivites_cache_file = os.path.join(self._cache_dir, f"{target_chembl_id}_aids.json")

        # 2) Fetch all activity records for this target, or load from cache
        self._logger.info(f"Fetching/loading all activities for ChEMBL ID {target_chembl_id}...")
        all_activity_records = self._get_cached_or_fetch(
            cache_file_path=pubchem_acitivites_cache_file,
            fetch_function=lambda: self._chembl_api_client.get_activities_for_target(
                target_chembl_id, self._bioactivity_measures
            ),
            data_type='ChEMBL activity records',
            force_refresh=force_refresh
        )
        self._logger.info(f"Found {len(all_activity_records)} total activity records.")

        # 3) Group all activity records by compound ID
        grouped_by_compound = defaultdict(list)
        for record in all_activity_records:
            chembl_id = record.get('molecule_chembl_id')
            if chembl_id:
                grouped_by_compound[chembl_id].append(record)

        # 4) Process each unique compound to calculate stats and create final object
        all_bioactives: List[BioactiveCompound] = []
        for chembl_id, records in grouped_by_compound.items():

            # 4.1) Group this compound's activities by measure type, converting units to nM
            grouped_activities = defaultdict(list)
            for record in records:
                unit = str(record.get('standard_units', '')).upper()
                value = record.get('standard_value')
                activity_type = str(record.get('standard_type', '')).upper()

                if not value:
                    continue

                conversion_factor = CONVERSION_FACTORS_TO_NM.get(unit)
                if conversion_factor:
                    try:
                        value_nm = float(value) * conversion_factor
                        grouped_activities[activity_type].append(value_nm)
                    except (ValueError, TypeError):
                        continue

            final_measure_type = None
            final_values = []
            for measure in self._bioactivity_measures:
                measure_upper = measure.upper()
                if grouped_activities[measure_upper]:
                    final_measure_type = measure_upper
                    final_values = grouped_activities[measure_upper]
                    break

            if not final_values:
                continue

            # 4.2) Calculate bioassay data statistics
            count = len(final_values)
            compound_bioassay_data = {
                "activity_type": final_measure_type,
                "activity_value": min(final_values),
                "n_measurements": count,
                "mean_activity": statistics.mean(final_values) if count > 0 else None,
                "median_activity": statistics.median(final_values) if count > 0 else None,
                "std_dev_activity": statistics.stdev(final_values) if count > 1 else 0.0,
            }

            # 4.3) Create the final BioactiveCompound object using data from the first record
            #     (since molecule properties will be the same across all records for this compound)
            first_record = records[0]
            canonical_smiles = first_record.get('canonical_smiles')
            data_validity_comment = first_record.get('data_validity_comment', None)
            inchikey = None
            mol_formula = None
            mol_weight = None

            if (data_validity_comment and data_validity_comment.upper() == CHEMBL_INVALID_DATA_COMMENT) or not canonical_smiles:
                continue

            mol = Chem.MolFromSmiles(canonical_smiles)
            if mol is not None:
                inchikey = Chem.MolToInchiKey(mol)
                mol_formula = CalcMolFormula(mol)
                mol_weight = Descriptors.MolWt(mol)

            format_mol_weight = None
            if mol_weight is not None:
                try:
                    format_mol_weight = round(float(mol_weight), 2)
                except (ValueError, TypeError):
                    format_mol_weight = None

            compound_obj = BioactiveCompound(
                source_db="ChEMBL",
                source_id=chembl_id,
                smiles=canonical_smiles,
                source_inchikey=inchikey,
                iupac_name=first_record.get('iupac_name', None),
                molecular_formula=mol_formula,
                molecular_weight=format_mol_weight,
                raw_data=records,  # Store all records for this compound
                **compound_bioassay_data  # Unpack the statistics dictionary
            )
            all_bioactives.append(compound_obj)

        # 5) Filter the final list by the 'activity_value' if a threshold was provided.
        if self._bioactivity_threshold is not None:
            self._logger.info(
                f"Filtering {len(all_bioactives)} compounds with threshold: <= {self._bioactivity_threshold} nM"
            )
            filtered_bioactives = [
                compound for compound in all_bioactives if compound.activity_value <= self._bioactivity_threshold
            ]
            self._logger.info(f"Found {len(filtered_bioactives)} compounds after filtering.")

            return filtered_bioactives

        return all_bioactives


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
        cache_dir: str = './data/cache',
        bioactivity_threshold: Optional[float] = None,
        # In nM (e.g. 1000 nM threshold to filter for compounds with Kd <= 1 µM)
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(bioactivity_measures, bioactivity_threshold, cache_dir, logger)
        self._api_client = PubChemAPIClient(logger=self._logger)

    def get_bioactive_compounds(self, target_uniprot_id: str, force_refresh: bool = False, ) -> List[BioactiveCompound]:
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

        # 2) Fetch AID list using the using API, or load from cache
        os.makedirs(self._cache_dir, exist_ok=True)
        aids_cache_file = os.path.join(self._cache_dir, f"{target_gene_id}_aids.json")
        aid_list = self._get_cached_or_fetch(
            cache_file_path=aids_cache_file,
            fetch_function=lambda: self._api_client.get_active_aids(target_gene_id),
            data_type='PubChem AIDs',
            force_refresh=force_refresh
        )

        if not aid_list:
            self._logger.warning(f"No assay IDs (AIDs) found for GeneID {target_gene_id}.")
            return []

        # 3) For each assay, retrieve active compound IDs (CIDs) and aggregate them.
        cids_cache_file = os.path.join(self._cache_dir, f"{target_gene_id}_cids.json")
        active_cids_list = self._get_cached_or_fetch(
            cache_file_path=cids_cache_file,
            fetch_function=lambda: self._fetch_all_cids(aids_list=aid_list),
            data_type='PubChem CIDs',
            force_refresh=force_refresh
        )

        if not active_cids_list:
            self._logger.error(f"No active compounds found for GeneID {target_gene_id}.")
            return []

        # 4) Retrieve full `pubchempy.Compound` objects for all CIDs.
        pubchempy_compound_api_start: float = time.time()
        pubchempy_compound_cache_file = os.path.join(self._cache_dir, f"{target_gene_id}_pubchempy_compounds.pkl")
        pubchempy_compounds = self._get_cached_or_fetch(
            cache_file_path=pubchempy_compound_cache_file,
            fetch_function=lambda: get_compounds_in_batches(cids=active_cids_list, logger=self._logger),
            data_type='PubChem bioactive `pubchempy` compound',
            use_pickle=True,
            force_refresh=force_refresh
        )
        pubchempy_compound_api_end: float = time.time()
        self._logger.info(f'PubChem bioactive compounds from CIDs total API query time: '
                          f'{round(pubchempy_compound_api_end - pubchempy_compound_api_start)} seconds')

        # 5) Fetch bioassay data for all retrieved compounds.
        bioassay_cache_file = os.path.join(self._cache_dir, f"{target_gene_id}_cid_bioassay_map.json")
        cid_to_bioassay_map = self._get_cached_or_fetch(
            cache_file_path=bioassay_cache_file,
            fetch_function=lambda: self._fetch_all_compound_bioassays(
                pubchempy_compounds=pubchempy_compounds,
                target_gene_id=target_gene_id
            ),
            data_type='PubChem Compound Bioassay',
            force_refresh=force_refresh
        )

        # 6) Create complete list of `BioactiveCompound` objects
        bioactivecompound_cache_file = os.path.join(
            self._cache_dir,
            f"{target_gene_id}_unfiltered_bioactivecompounds.pkl"
        )
        all_bioactives: List[BioactiveCompound] = self._get_cached_or_fetch(
            cache_file_path=bioactivecompound_cache_file,
            fetch_function=lambda: self._get_all_bioactive_compounds(
                pubchempy_compounds=pubchempy_compounds,
                cid_to_bioassay_map=cid_to_bioassay_map
            ),
            data_type='BioactiveCompound object',
            use_pickle=True,
            force_refresh=force_refresh
        )

        # 7) Filter final list of `BioactiveCompound` objects by potency if threshold is provided.
        if self._bioactivity_threshold is not None:
            self._logger.info(f"Filtering {len(all_bioactives)} compounds with threshold: "
                              f"<= {self._bioactivity_threshold} nM")
            filtered_bioactives: List[BioactiveCompound] = [
                compound for compound in all_bioactives if compound.activity_value <= self._bioactivity_threshold
            ]
            self._logger.info(f"Found {len(filtered_bioactives)} compounds after filtering.")

            return filtered_bioactives

        return all_bioactives

    def _fetch_all_cids(self, aids_list: List[int]) -> List[int]:
        cids_api_start: float = time.time()
        active_cids = set()

        # Create a new partial function with `logger` argument fixed. This allows us to pass a fixed `logger` argument
        # to the `get_active_cids_wrapper()` function when it is mapped to each AID element in `aid_list` via
        # `concurrent.futures.ThreadPoolExecutor.map()`
        get_active_cids_partial = partial(self._api_client.get_active_cids)

        # Create thread pool using Python’s `ThreadPoolExecutor` to issue multiple API calls concurrently in batches
        with concurrent.futures.ThreadPoolExecutor(max_workers=9) as executor:
            # Map and apply partial function of `cids_for_aid_wrapper()` to every element in `aid_list` concurrently
            results = list(executor.map(get_active_cids_partial, aids_list))

            for cids in results:
                active_cids.update(cids)

        cids_api_end: float = time.time()
        self._logger.info(f'PubChem CID total API query time: {round(cids_api_end - cids_api_start)} seconds')

        return list(active_cids)

    def _fetch_all_compound_bioassays(
        self,
        pubchempy_compounds: List[pcp.Compound],
        target_gene_id: str
    ) -> Dict[str, Any]:
        self._logger.info(f"Fetching bioassay data for {len(pubchempy_compounds)} compounds...")
        potencies_api_start: float = time.time()
        cid_to_bioassay_map = {}

        # Create a new partial function with `target_gene_id` and `logger` argument fixed. As before, this allows
        # us to pass these fixed arguments to `self._get_compound_bioassay_data()` when it is mapped to each
        # compound element in the batched `bioactive_compounds` iterable via
        # `concurrent.futures.ThreadPoolExecutor.map()`
        get_compound_bioassay_data_partial = partial(
            self._api_client.get_compound_bioassay_data,
            target_gene_id=target_gene_id,
            bioactivity_measures=self._bioactivity_measures
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
                        cid_to_bioassay_map[compound.cid] = bioassay_data

        potencies_api_end: float = time.time()
        self._logger.info(f"PubChem bioactive compound bioassays total API query time: "
                          f"{round(potencies_api_end - potencies_api_start)} seconds\n"
                          f"Found bioassay data for {len(cid_to_bioassay_map)} compounds.")

        return cid_to_bioassay_map

    def _get_all_bioactive_compounds(
        self,
        pubchempy_compounds: List[pcp.Compound],
        cid_to_bioassay_map: Dict[str, Any]
    ) -> List[BioactiveCompound]:
        all_bioactives: List[BioactiveCompound] = []
        for pubchempy_compound in pubchempy_compounds:
            compound_bioassay_data: Dict = cid_to_bioassay_map.get(str(pubchempy_compound.cid))

            if compound_bioassay_data is None:
                self._logger.debug(f"Skipping compound CID {pubchempy_compound.cid} due to missing bioassay data.")
                continue

            compound_obj = self._create_bioactive_compound(
                pubchempy_compound=pubchempy_compound,
                bioassay_data=compound_bioassay_data
            )
            if compound_obj:
                all_bioactives.append(compound_obj)

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
        smiles = getattr(pubchempy_compound, 'canonical_smiles', None)
        if not smiles:
            return None

        mol = Chem.MolFromSmiles(smiles)
        source_inchikey = getattr(pubchempy_compound, 'inchikey', None)
        source_mol_formula = getattr(pubchempy_compound, 'molecular_formula', None)
        source_mol_weight = getattr(pubchempy_compound, 'molecular_weight', None)

        if mol:
            final_inchikey = source_inchikey or Chem.MolToInchiKey(mol)
            final_mol_formula = source_mol_formula or CalcMolFormula(mol)
            final_mol_weight = source_mol_weight or Descriptors.MolWt(mol)
        else:
            final_inchikey = source_inchikey
            final_mol_formula = source_mol_formula
            final_mol_weight = source_mol_weight

        format_mol_weight = None
        if final_mol_weight is not None:
            try:
                format_mol_weight = round(float(final_mol_weight), 2)
            except (ValueError, TypeError):
                format_mol_weight = None

        return BioactiveCompound(
            source_db='PubChem',
            source_id=str(pubchempy_compound.cid),
            smiles=smiles,
            activity_type=bioassay_data['activity_type'],
            activity_value=bioassay_data['best_value'],
            source_inchikey=final_inchikey,
            iupac_name=getattr(pubchempy_compound, 'iupac_name', None),
            molecular_formula=final_mol_formula,
            molecular_weight=format_mol_weight,
            n_measurements=bioassay_data["n_measurements"],
            mean_activity=round(bioassay_data["mean_value"], 2) if bioassay_data["mean_value"] else None,
            median_activity=round(bioassay_data["median_value"], 2) if bioassay_data["median_value"] else None,
            std_dev_activity=round(bioassay_data["std_dev_value"], 2) if bioassay_data["std_dev_value"] else 0.0,
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
