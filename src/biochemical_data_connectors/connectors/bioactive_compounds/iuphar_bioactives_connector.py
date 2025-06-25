import os
import statistics
import logging
from collections import defaultdict
from typing import Optional, List

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from biochemical_data_connectors.connectors.bioactive_compounds.base_bioactives_connector import BaseBioactivesConnector
from biochemical_data_connectors.utils.api.iuphar_api import IupharApiClient
from biochemical_data_connectors.models import BioactiveCompound
from biochemical_data_connectors.utils.files_utils import get_cached_or_fetch
from biochemical_data_connectors.utils.standardization_utils import convert_p_value_to_nm


class IUPHARBioactivesConnector(BaseBioactivesConnector):
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,
        cache_dir: str = '/data/cache',
        logger: Optional[logging.Logger] = None
    ):
        super().__init__(
            bioactivity_measures=bioactivity_measures,
            bioactivity_threshold=bioactivity_threshold,
            cache_dir=cache_dir,
            logger=logger
        )
        self._iuphar_api_client = IupharApiClient(logger=logger)

    def get_bioactive_compounds(
        self,
        target_uniprot_id: str,
        force_refresh: bool = False
    ) -> List[BioactiveCompound]:
        # 1. Query for the target ID by UniProt ID.
        iuphar_target_id: int = self._iuphar_api_client.get_iuphar_target_id(uniprot_id=target_uniprot_id)
        if iuphar_target_id is None:
            self._logger.error(f"No matching IUPHAR-BPS target found for UniProt ID {target_uniprot_id}")
            return []

        # 2. Fetch all activity records for this target, using the cache if available.
        os.makedirs(self._cache_dir, exist_ok=True)
        iuphar_activities_cache_file = os.path.join(self._cache_dir, f"IUPHAR-BPS/{iuphar_target_id}.json")

        self._logger.info(f"Fetching/loading all IUPHAR-BPS activities for Uniprot ID {target_uniprot_id}...")
        p_activity_measures =  ['p' + measure for measure in self._bioactivity_measures]
        all_iuphar_activity_records = get_cached_or_fetch(
            cache_file_path=iuphar_activities_cache_file,
            fetch_function=lambda: self._iuphar_api_client.get_actives_from_target_id(
                target_id=iuphar_target_id,
                p_bioactivity_measures=p_activity_measures,
            ),
            data_type='IUPHAR-BPS activity records',
            force_refresh=force_refresh,
            logger=self._logger
        )

        # 3. Group all activity records by ligand ID
        grouped_by_compound = defaultdict(list)
        for record in all_iuphar_activity_records:
            ligand_id = record.get('ligandId')
            if ligand_id:
                grouped_by_compound[ligand_id].append(record)

        # 4. Process each unique compound to calculate stats and create final object
        all_bioactives: List[BioactiveCompound] = []
        for ligand_id, records in grouped_by_compound.items():
            first_record = records[0]

            # 4.1. Collect all valid, numeric affinity values for this compound.
            #      The BindingDB API returns values in nM, so no conversion is needed.
            final_values = []
            p_measure = first_record.get('affinityParameter')
            final_measure_type = p_measure.replace(p_measure[0], '', 1)
            for record in records:
                try:
                    value = round(convert_p_value_to_nm(float(record.get('affinity'))), 2)
                    final_values.append(value)
                except (ValueError, TypeError):
                    continue

            if not final_values:
                continue

            # 4.2. Calculate bioassay data statistics
            count = len(final_values)
            stats = {
                "activity_type": final_measure_type,
                "activity_value": min(final_values),
                "n_measurements": count,
                "mean_activity": round(statistics.mean(final_values), 2) if count > 0 else None,
                "median_activity": round(statistics.median(final_values), 2) if count > 0 else None,
                "std_dev_activity": round(statistics.stdev(final_values), 2) if count > 1 else 0.0,
            }

            # 4.3. BindingDB response doesn't provide InCHIKey, molecular formula, or molecular weight.
            #      Use RDKit to calculate these for consistency and create final BioactiveCompound object.
            #      N.B. The API call already filtered by threshold, so we don't need to filter again.
            mol_data = self._iuphar_api_client.get_mol_data_from_ligand_id(ligand_id=ligand_id)
            if not mol_data.get('smiles'):
                continue

            mol = Chem.MolFromSmiles(mol_data.get('smiles'))
            if not mol:
                continue

            compound_obj = BioactiveCompound(
                source_db="IUPHAR/BPS Guide to PHARMACOLOGY",
                source_id=ligand_id,
                smiles=mol_data.get('smiles'),
                target_uniprot=target_uniprot_id,
                source_inchikey=mol_data.get('inchikey') if mol_data.get('inchikey') else Chem.MolToInchiKey(mol),
                molecular_formula=CalcMolFormula(mol),
                molecular_weight=mol_data.get('molecular_weight') if mol_data.get('molecular_weight') else round(
                    Descriptors.MolWt(mol), 2
                ),
                raw_data=records,
                **stats
            )
            all_bioactives.append(compound_obj)

        if self._bioactivity_threshold is not None:
            self._logger.info(
                f"Filtering {len(all_bioactives)} IUPHAR/BPS compounds with threshold: <= {self._bioactivity_threshold} nM"
            )
            filtered_bioactives = [
                compound for compound in all_bioactives if compound.activity_value <= self._bioactivity_threshold
            ]
            self._logger.info(f"Found {len(filtered_bioactives)} IUPHAR/BPS compounds after filtering.")

            return filtered_bioactives

        return all_bioactives
