import logging
import time
import requests
from typing import Dict, Any, Optional, List

from biochemical_data_connectors.utils.api.base_api import BaseApiClient
from biochemical_data_connectors.constants import RestApiEndpoints


class IupharApiClient(BaseApiClient):
    def __init__(self, logger: Optional[logging.Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_iuphar_target_id(self, uniprot_id: str) -> Optional[int]:
        """
        Queries the IUPHAR/BPS API for the internal target ID.

        Maps a UniProt accession ID to the corresponding internal database ID
        used by the Guide to PHARMACOLOGY.

        Parameters
        ----------
        uniprot_id : str
            The UniProt accession ID of the target.

        Returns
        -------
        Optional[int]
            The internal integer target ID if found, otherwise None.
        """
        iuphar_target_query_start = time.time()
        iuphar_target_url = RestApiEndpoints.IUPHAR_TARGET_ID_FROM_UNIPROT.url(
            uniprot_id=uniprot_id
        )
        self._logger.info(f'Querying IUPHAR/BPS Guide to Pharmacology API for Uniprot {uniprot_id} target ID')
        try:
            response = self._session.get(iuphar_target_url, timeout=15)
            response.raise_for_status()
            target_id = response.json()[0].get('targetId')

            iuphar_target_query_end = time.time()
            self._logger.info(
                f'IUPHAR total query time: {round(iuphar_target_query_end - iuphar_target_query_start)} seconds'
            )
            if not target_id:
                return None

            return target_id

        except requests.exceptions.RequestException as e:
            self._logger.error(f'IUPHAR Guide to Pharmacology API request failed for target {uniprot_id}: {e}')

            return None

    def get_actives_from_target_id(
        self,
        target_id: int,
        p_bioactivity_measures: List[str],
    ) -> List[Any]:
        """
        Queries IUPHAR/BPS API for interactions and returns standardized data.

        This method fetches interactions for a target, filtering by the desired
        p-value activity types at the API level.

        Parameters
        ----------
        target_id : int
            The internal IUPHAR/BPS target ID.
        p_bioactivity_measures : List[str]
            A list of p-value activity types to fetch (e.g., ['pKi', 'pIC50']).

        Returns
        -------
        List[Dict]
            A list of dictionaries, each representing a clean interaction record.
        """
        all_records = []
        iuphar_interactions_query_start = time.time()

        # 1. Iterate through the desired measure types and make a separate API call for each.
        for p_measure in p_bioactivity_measures:
            self._logger.info(f'Querying IUPHAR/BPS for {p_measure} data for target ID {target_id}...')

            # 2. Build the URL with the affinityType filter.
            iuphar_interactions_url = RestApiEndpoints.IUPHAR_INTERACTIONS_FROM_TARGET_ID_FILTERED.url(
                target_id=target_id,
                p_measure=p_measure
            )
            try:
                response = self._session.get(iuphar_interactions_url, timeout=15)
                response.raise_for_status()
                iuphar_interactions = response.json()

                # 3. Defensive logic to ensure interactions are returned and contain affinity
                if not iuphar_interactions:
                    self._logger.warning(f'No active ligand interactions found for target ID {target_id}')
                    continue

                for interaction in iuphar_interactions:
                    if interaction.get('affinity') is None:
                        continue

                    all_records.append(interaction)
            except requests.exceptions.RequestException as e:
                self._logger.error(f'Error querying IUPHAR/BPS for {p_measure} at target {target_id}: {e}')

        iuphar_interactions_query_end = time.time()
        self._logger.info(
            f'IUPHAR interactions total query time: '
            f'{round(iuphar_interactions_query_end - iuphar_interactions_query_start)} seconds'
        )
        self._logger.info(f'IUPHAR/BPS query found a total of {len(all_records)} relevant records.')

        return all_records


    def get_mol_data_from_ligand_id(self, ligand_id: str) -> Dict:
        """
        Fetches structural and property data for a given IUPHAR ligand ID.

        This method makes two separate API calls to retrieve structure (SMILES, etc.)
        and molecular properties (molecular weight) and combines them.

        Parameters
        ----------
        ligand_id : str
            The internal IUPHAR/BPS ligand ID.

        Returns
        -------
        Dict[str, Any]
            A dictionary containing the combined molecular data.
        """
        iuphar_ligand_structure_url = RestApiEndpoints.IUPHAR_LIGAND_STRUCTURE_FROM_LIGAND_ID.url(
            ligand_id=ligand_id
        )
        iuphar_ligand_mol_props_url = RestApiEndpoints.IUPHAR_LIGAND_MOLECULAR_PROPERTIES_FROM_LIGAND_ID.url(
            ligand_id=ligand_id
        )
        mol_data: Dict = {}
        try:
            response = self._session.get(iuphar_ligand_structure_url, timeout=15)
            response.raise_for_status()
            iuphar_ligand_structure = response.json()

            if not iuphar_ligand_structure:
                self._logger.warning(f'No structural data found for ligand ID {ligand_id}')
                mol_data |= {'smiles': None, 'inchikey': None, 'iupac_name': None}

            mol_data |= {
                'smiles': iuphar_ligand_structure.get('smiles'),
                'inchikey': iuphar_ligand_structure.get('inchiKey'),
                'iupac_name': iuphar_ligand_structure.get('iupacName')
            }
        except requests.exceptions.RequestException as e:
            self._logger.error(f"Error querying IUPHAR/BPS ligand {ligand_id} structure: {e}")
            mol_data |= {'smiles': None, 'inchikey': None, 'iupac_name': None}

        try:
            response = self._session.get(iuphar_ligand_mol_props_url, timeout=15)
            response.raise_for_status()
            iuphar_mol_properties = response.json()

            if not iuphar_mol_properties or not iuphar_mol_properties.get('molecularWeight'):
                self._logger.warning(f'No molecular property data found for ligand ID {ligand_id}')
                mol_data |= {'molecular_weight': None}

            mol_data |= {'molecular_weight': round(iuphar_mol_properties.get('molecularWeight'), 2)}
        except requests.exceptions.RequestException as e:
            self._logger.error(f'Error querying IUPHAR/BPS ligand {ligand_id} molecular properties: {e}')
            mol_data |= {'molecular_weight': None}

        return mol_data
