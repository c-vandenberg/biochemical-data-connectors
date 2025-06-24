import time
import requests
import logging
from typing import List, Optional

from biochemical_data_connectors.utils.api.base_api import BaseAPIClient
from biochemical_data_connectors.constants import RestApiEndpoints


class BindingDBAPIClient(BaseAPIClient):
    def __init__(self, logger: Optional[logging.Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_actives_from_target(
        self,
        uniprot_id: str,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,  # In nM.
    ) -> List:
        bindingdb_start = time.time()

        cutoff_str = f";{int(bioactivity_threshold)}" if bioactivity_threshold is not None else ""
        url = RestApiEndpoints.BINDINGDB_LIGANDS_FROM_UNIPROT_ID.url(
            uniprot_id=uniprot_id,
            cutoff_str=cutoff_str
        )
        self._logger.info(f"Querying BindingDB for target: {uniprot_id}")
        try:
            response = self._session.get(url)
            response.raise_for_status()
            data = response.json()
            bdb_affinities = data.get('getLindsByUniprotResponse', {}).get('bdb.affinities', [])
            if not bdb_affinities:
                self._logger.warning(f"No BindingDB actives found for {uniprot_id}.")
                return []

            filtered_bdb_actives = [
                record for record in bdb_affinities if record.get('bdb.affinity_type') in bioactivity_measures
            ]
            bindingdb_end = time.time()
            self._logger.info(f'Binding DB total query time: {round(bindingdb_end - bindingdb_start)} seconds')

            return filtered_bdb_actives

        except requests.exceptions.RequestException as e:
            self._logger.error(f"Error querying BindingDB: {e}")

            return []