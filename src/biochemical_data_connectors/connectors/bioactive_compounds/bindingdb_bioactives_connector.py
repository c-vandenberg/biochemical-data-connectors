import os
import statistics
import logging
from collections import defaultdict
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from biochemical_data_connectors.utils.api.bindingdb_api import BindingDbApiClient
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.connectors.bioactive_compounds.base_bioactives_connector import BaseBioactivesConnector
from biochemical_data_connectors.utils.files_utils import get_cached_or_fetch


class BindingDbBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts and processes bioactive compounds from BindingDB.

    This connector orchestrates the fetching of data for a given UniProt ID,
    groups all measurements for each unique compound, calculates statistics,
    and returns a standardized list of `BioactiveCompound` objects.

    Attributes
    ----------
    _bdb_api_client : BindingDbApiClient
        An instance of the client used to handle all direct API communications.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,  # In nM.
        cache_dir: str = './data/cache',
        logger: Optional[logging.Logger] = None,
    ):
        super().__init__(
            bioactivity_measures=bioactivity_measures,
            bioactivity_threshold=bioactivity_threshold,
            cache_dir=cache_dir,
            logger=logger
        )
        self._bdb_api_client = BindingDbApiClient(logger=self._logger)

    def get_bioactive_compounds(
        self,
        target_uniprot_id: str,
        force_refresh: bool = False
    ) -> List[BioactiveCompound]:
        """
        Retrieve and process bioactive compound data for a target from BindingDB.

        This method orchestrates the entire workflow:
        1. Fetches all relevant activity records from the BindingDB API, using a
           cache to improve performance on subsequent runs.
        2. Groups all records by their unique compound `bdb.monomerid`.
        3. For each unique compound, calculates statistics on its affinity values.
        4. Creates a standardized `BioactiveCompound` object containing the
           aggregated data and RDKit-calculated properties.

        Parameters
        ----------
        target_uniprot_id : str
            The UniProt accession for the target.
        force_refresh : bool, optional
            If True, ignores any existing cache and forces a new API call.

        Returns
        -------
        List[BioactiveCompound]
            A list of fully populated and standardized BioactiveCompound objects.
        """
        # 1. Fetch all activity records for this target, using the cache if available.
        os.makedirs(self._cache_dir, exist_ok=True)
        bdb_activities_cache_file = os.path.join(self._cache_dir, f'BindingDB/{target_uniprot_id}.json')

        self._logger.info(f'Fetching/loading all BindingDB activities for Uniprot ID {target_uniprot_id}...')
        all_bdb_activity_records = get_cached_or_fetch(
            cache_file_path=bdb_activities_cache_file,
        fetch_function=lambda: self._bdb_api_client.get_actives_from_target_uniprot(
                uniprot_id=target_uniprot_id,
                bioactivity_measures=self._bioactivity_measures,
                bioactivity_threshold=self._bioactivity_threshold
            ),
            data_type='BindingDB activity records',
            force_refresh=force_refresh,
            logger=self._logger
        )

        if not all_bdb_activity_records:
            return []

        # 2. Group all activity records by BindingDB monomer ID
        grouped_by_compound = defaultdict(list)
        for record in all_bdb_activity_records:
            monomer_id = record.get('bdb.monomerid')
            if monomer_id:
                grouped_by_compound[monomer_id].append(record)

        # 3. Process each unique compound to calculate stats and create a final object.
        all_bioactives: List[BioactiveCompound] = []
        for monomer_id, records in grouped_by_compound.items():
            first_record = records[0]

            # 3.1. Collect all valid, numeric affinity values for this compound.
            #      The BindingDB API returns values in nM, so no conversion is needed.
            final_values = []
            final_measure_type = first_record.get('bdb.affinity_type')
            for record in records:
                try:
                    value = float(record.get('bdb.affinity'))
                    final_values.append(value)
                except (ValueError, TypeError):
                    continue

            if not final_values:
                continue

            # 3.2. Calculate bioassay data statistics
            count = len(final_values)
            stats = {
                'activity_type': final_measure_type,
                'activity_value': min(final_values),
                'n_measurements': count,
                'mean_activity': round(statistics.mean(final_values), 2) if count > 0 else None,
                'median_activity': round(statistics.median(final_values), 2) if count > 0 else None,
                'std_dev_activity': round(statistics.stdev(final_values), 2) if count > 1 else 0.0,
            }

            # 3.3. BindingDB response doesn't provide InCHIKey, molecular formula, or molecular weight.
            #      Use RDKit to calculate these for consistency and create final BioactiveCompound object.
            #      N.B. The API call already filtered by threshold, so we don't need to filter again.
            mol = Chem.MolFromSmiles(first_record.get('bdb.smile'))
            if not mol:
                continue

            compound_obj = BioactiveCompound(
                source_db='BindingDB',
                source_id=monomer_id,
                smiles=Chem.MolToSmiles(mol, canonical=True),
                target_uniprot=target_uniprot_id,
                source_inchikey=Chem.MolToInchiKey(mol),
                molecular_formula=CalcMolFormula(mol),
                molecular_weight=round(Descriptors.MolWt(mol), 2),
                raw_data=records,
                **stats
            )
            all_bioactives.append(compound_obj)

        return all_bioactives
