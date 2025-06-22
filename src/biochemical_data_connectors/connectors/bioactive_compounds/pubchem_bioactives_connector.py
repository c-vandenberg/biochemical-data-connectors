import os
import time
import logging
import concurrent.futures
from functools import partial
from typing import List, Dict, Optional, Any

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from biochemical_data_connectors.connectors.bioactive_compounds.base_bioactives_connector import BaseBioactivesConnector
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.utils.iter_utils import batch_iterable
from biochemical_data_connectors.utils.api.pubchem_api import PubChemAPIClient, get_compounds_in_batches
from biochemical_data_connectors.utils.api.mappings import uniprot_to_gene_id_mapping
from biochemical_data_connectors.utils.files_utils import get_cached_or_fetch

CHEMBL_INVALID_DATA_COMMENT = 'OUTSIDE TYPICAL RANGE'


class PubChemBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts bioactive compounds for a given target from PubChem using a UniProt accession.

    This connector orchestrates a multi-step, cached workflow to query PubChem,
    retrieve all relevant compounds and their bioactivity data, and format
    them into standardized `BioactiveCompound` objects.

    Attributes
    ----------
    _api_client : PubChemAPIClient
        An instance of the client used to handle all direct API communications.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None, # In nM.
        cache_dir: str = './data/cache',
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(
            bioactivity_measures=bioactivity_measures,
            bioactivity_threshold=bioactivity_threshold,
            cache_dir=cache_dir,
            logger=logger
        )
        self._api_client = PubChemAPIClient(logger=self._logger)

    def get_bioactive_compounds(
        self,
        target_uniprot_id: str,
        force_refresh: bool = False
    ) -> List[BioactiveCompound]:
        """
        Retrieve and process bioactive compound data for a target from PubChem.
        The target is provided as a UniProt accession (e.g. "P00533").

        This method orchestrates the entire data fetching pipeline, using a
        multi-level cache to accelerate subsequent runs. Each major step
        (fetching AIDs, CIDs, Compound objects, etc.) is cached independently.

        This pipeline consists of the following steps:

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
        force_refresh : bool, optional
            If True, ignores all caches and re-fetches all data from the APIs.

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

        # 2) Fetch Assay ID (AID) list using PubChem API, or load from cache
        os.makedirs(self._cache_dir, exist_ok=True)
        aids_cache_file = os.path.join(self._cache_dir, f"pubchem/{target_gene_id}_aids.json")
        aid_list = get_cached_or_fetch(
            cache_file_path=aids_cache_file,
            fetch_function=lambda: self._api_client.get_active_aids(target_gene_id),
            data_type='PubChem AIDs',
            force_refresh=force_refresh,
            logger=self._logger
        )

        if not aid_list:
            self._logger.warning(f"No assay IDs (AIDs) found for GeneID {target_gene_id}.")
            return []

        # 3) Fetch active compound IDs (CIDs) for each assay using PubChem API, or load from cache.
        cids_cache_file = os.path.join(self._cache_dir, f"pubchem/{target_gene_id}_cids.json")
        active_cids_list = get_cached_or_fetch(
            cache_file_path=cids_cache_file,
            fetch_function=lambda: self._fetch_all_cids(aids_list=aid_list),
            data_type='PubChem CIDs',
            force_refresh=force_refresh,
            logger=self._logger
        )

        if not active_cids_list:
            self._logger.error(f"No active compounds found for GeneID {target_gene_id}.")
            return []

        # 4) Fetch full `pubchempy.Compound` objects for all CIDs using PubChem API, or load from cache.
        pubchempy_compound_api_start: float = time.time()
        pubchempy_compound_cache_file = os.path.join(
            self._cache_dir,
            f"pubchem/{target_gene_id}_pubchempy_compounds.pkl"
        )
        pubchempy_compounds = get_cached_or_fetch(
            cache_file_path=pubchempy_compound_cache_file,
            fetch_function=lambda: get_compounds_in_batches(cids=active_cids_list, logger=self._logger),
            data_type='PubChem bioactive `pubchempy` compound',
            use_pickle=True,
            force_refresh=force_refresh,
            logger=self._logger
        )
        pubchempy_compound_api_end: float = time.time()
        self._logger.info(f'PubChem bioactive compounds from CIDs total API query time: '
                          f'{round(pubchempy_compound_api_end - pubchempy_compound_api_start)} seconds')

        # 5) Fetch bioassay data for all `pubchempy.Compound` compounds using PubChem API, or load from cache.
        bioassay_cache_file = os.path.join(
            self._cache_dir,
            f"pubchem/{target_gene_id}_cid_bioassay_map.json"
        )
        cid_to_bioassay_map = get_cached_or_fetch(
            cache_file_path=bioassay_cache_file,
            fetch_function=lambda: self._fetch_all_compound_bioassays(
                pubchempy_compounds=pubchempy_compounds,
                target_gene_id=target_gene_id
            ),
            data_type='PubChem Compound Bioassay',
            force_refresh=force_refresh,
            logger=self._logger
        )

        # 6) Create the final list of `BioactiveCompound` objects, using cache if available.
        bioactivecompound_cache_file = os.path.join(
            self._cache_dir,
            f"pubchem/{target_gene_id}_unfiltered_bioactivecompounds.pkl"
        )
        all_bioactives: List[BioactiveCompound] = get_cached_or_fetch(
            cache_file_path=bioactivecompound_cache_file,
            fetch_function=lambda: self._get_all_bioactive_compounds(
                target_uniprot_id=target_uniprot_id,
                pubchempy_compounds=pubchempy_compounds,
                cid_to_bioassay_map=cid_to_bioassay_map
            ),
            data_type='BioactiveCompound object',
            use_pickle=True,
            force_refresh=force_refresh,
            logger=self._logger
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
        """
        Concurrently fetches all active CIDs for a given list of AIDs.

        Parameters
        ----------
        aids_list : List[str]
            A list of PubChem Assay IDs (AIDs) to process.

        Returns
        -------
        List[int]
            A deduplicated list of active Compound IDs (CIDs).
        """
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
        """
        Concurrently fetches bioassay data for a list of compounds.

        Parameters
        ----------
        pubchempy_compounds : List[pcp.Compound]
            The list of compounds to fetch data for.
        target_gene_id : str
            The NCBI GeneID to filter assays by.

        Returns
        -------
        Dict[int, Dict[str, Any]]
            A dictionary mapping a CID to its fetched bioassay data.
        """
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
                        cid_to_bioassay_map[str(compound.cid)] = bioassay_data

        potencies_api_end: float = time.time()
        self._logger.info(f"PubChem bioactive compound bioassays total API query time: "
                          f"{round(potencies_api_end - potencies_api_start)} seconds\n"
                          f"Found bioassay data for {len(cid_to_bioassay_map)} compounds.")

        return cid_to_bioassay_map

    def _get_all_bioactive_compounds(
        self,
        target_uniprot_id: str,
        pubchempy_compounds: List[pcp.Compound],
        cid_to_bioassay_map: Dict[str, Any]
    ) -> List[BioactiveCompound]:
        """
        Assembles the final list of `BioactiveCompound` objects.

        This method iterates through the fetched compound data, combines it with
        the corresponding bioassay data, and creates the final,
        standardized `BioactiveCompound` objects.

        Parameters
        ----------
        pubchempy_compounds : List[pcp.Compound]
            The list of raw compound objects from PubChem.
        cid_to_bioassay_map : Dict[int, Dict[str, Any]]
            A mapping of CIDs to their fetched bioassay statistics.

        Returns
        -------
        List[BioactiveCompound]
            A list of the fully populated `BioactiveCompound` objects.
        """
        all_bioactives: List[BioactiveCompound] = []
        for pubchempy_compound in pubchempy_compounds:
            compound_bioassay_data: Dict = cid_to_bioassay_map.get(str(pubchempy_compound.cid))

            if compound_bioassay_data is None:
                self._logger.debug(f"Skipping compound CID {pubchempy_compound.cid} due to missing bioassay data.")
                continue

            compound_obj = self._create_bioactive_compound(
                target_uniprot_id=target_uniprot_id,
                pubchempy_compound=pubchempy_compound,
                bioassay_data=compound_bioassay_data
            )
            if compound_obj:
                all_bioactives.append(compound_obj)

        return all_bioactives

    @staticmethod
    def _create_bioactive_compound(
        target_uniprot_id: str,
        pubchempy_compound: pcp.Compound,
        bioassay_data: Dict[str, Any]
    ) -> Optional[BioactiveCompound]:
        """
        Helper to convert a `pubchempy.Compound` to a `BioactiveCompound`.

        This method safely extracts attributes from the source object and uses
        them to instantiate the standardized `BioactiveCompound` dataclass.

        Parameters
        ----------
        target_uniprot_id : str
            The UniProt accession for the target.
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
            target_uniprot=target_uniprot_id,
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