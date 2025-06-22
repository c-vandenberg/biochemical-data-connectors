import os
import logging
import statistics
from collections import defaultdict
from typing import List, Optional

from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from biochemical_data_connectors.connectors.bioactive_compounds.base_bioactives_connector import BaseBioactivesConnector
from biochemical_data_connectors.constants import CONVERSION_FACTORS_TO_NM
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.utils.api.chembl_api import ChEMBLAPIClient
from biochemical_data_connectors.utils.files_utils import get_cached_or_fetch

CHEMBL_INVALID_DATA_COMMENT = 'OUTSIDE TYPICAL RANGE'


class ChEMBLBioactivesConnector(BaseBioactivesConnector):
    """
    Extracts bioactive compounds from ChEMBL using a target's UniProt accession.

    This class orchestrates the fetching and processing of data from the ChEMBL
    database, handling UniProt ID lookups, API pagination, data aggregation,
    and unit conversion to produce a standardized list of `BioactiveCompound`
    objects.

    Attributes
    ----------
    _chembl_webresource_client : object
        A client for the high-level ChEMBL API, used for target lookups.
    _chembl_api_client : ChEMBLAPIClient
        A client for the low-level ChEMBL REST API, used for activity fetching.
    """

    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,  # In nM.
        cache_dir: str = './data/cache',
        logger: Optional[logging.Logger] = None,
        core_chembl_client=None,
    ):
        super().__init__(
            bioactivity_measures=bioactivity_measures,
            bioactivity_threshold=bioactivity_threshold,
            cache_dir=cache_dir,
            logger=logger
        )
        self._chembl_webresource_client = core_chembl_client if core_chembl_client else new_client
        self._chembl_api_client: ChEMBLAPIClient = ChEMBLAPIClient(logger=self._logger)

    def get_bioactive_compounds(
        self,
        target_uniprot_id: str,
        force_refresh: bool = False
    ) -> List[BioactiveCompound]:
        """
        Retrieve bioassay data for bioactive compounds from ChEMBL using a target's UniProt accession.

        This method queries the ChEMBL activity API, fetching full records
        for compounds that match the target and bioactivity criteria, and
        returns them as a list of structured BioactiveCompound objects.

        This method orchestrates the entire workflow:
        1. Converts the UniProt ID to a ChEMBL Target ID.
        2. Fetches all relevant activity records from the ChEMBL API, using a cache
           to improve performance on subsequent runs.
        3. Groups all records by their unique compound ID.
        4. For each compound, calculates statistics (count, mean, median, stdev)
           for the highest-priority bioactivity type found.
        5. Creates a standardized `BioactiveCompound` object containing the
           aggregated data.
        6. Applies an optional potency filter to the final list.

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
        # 1) Search for the target by UniProt ID and retrieve the first matching result
        target_results = self._chembl_webresource_client.target.filter(target_components__accession=target_uniprot_id)
        if not target_results:
            self._logger.error(f"No matching target found for UniProt ID {target_uniprot_id}")
            return []

        target_chembl_id = target_results[0]['target_chembl_id']

        # 2) Fetch all activity records for this target, using the cache if available.
        os.makedirs(self._cache_dir, exist_ok=True)
        chembl_acitivites_cache_file = os.path.join(self._cache_dir, f"chembl/{target_chembl_id}_aids.json")

        self._logger.info(f"Fetching/loading all activities for ChEMBL ID {target_chembl_id}...")
        all_activity_records = get_cached_or_fetch(
            cache_file_path=chembl_acitivites_cache_file,
            fetch_function=lambda: self._chembl_api_client.get_activities_for_target(
                target_chembl_id,
                self._bioactivity_measures
            ),
            data_type='ChEMBL activity records',
            force_refresh=force_refresh,
            logger=self._logger
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
            # 4.1) Find the first available canonical SMILES from all records for this compound.
            #      If none found, skip this compound.
            canonical_smiles = next((r.get('canonical_smiles') for r in records if r.get('canonical_smiles')), None)
            if not canonical_smiles:
                continue

            # 4.2) Group this compound's activities by measure type, converting units to nM
            grouped_activities = defaultdict(list)
            for record in records:
                data_validity_comment = record.get('data_validity_comment')

                # Skip record if data record is invalid
                if data_validity_comment and data_validity_comment.upper() == CHEMBL_INVALID_DATA_COMMENT:
                    continue

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

            # 4.3) Find the highest-priority activity data for this compound.
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

            # 4.4) Calculate bioassay data statistics
            count = len(final_values)
            compound_bioassay_data = {
                "activity_type": final_measure_type,
                "activity_value": min(final_values),
                "n_measurements": count,
                "mean_activity": statistics.mean(final_values) if count > 0 else None,
                "median_activity": statistics.median(final_values) if count > 0 else None,
                "std_dev_activity": statistics.stdev(final_values) if count > 1 else 0.0,
            }

            # 4.5) ChEMBL response doesn't provide InCHIKey, molecular formula, or molecular weight.
            #      Use RDKit to calculate these for consistency.
            inchikey = None
            mol_formula = None
            mol_weight = None

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

            # 4.6) Create the final BioactiveCompound object
            compound_obj = BioactiveCompound(
                source_db="ChEMBL",
                source_id=chembl_id,
                smiles=canonical_smiles,
                target_uniprot=target_uniprot_id,
                source_inchikey=inchikey,
                iupac_name=records[0].get('iupac_name', None),
                molecular_formula=mol_formula,
                molecular_weight=format_mol_weight,
                raw_data=records,
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
