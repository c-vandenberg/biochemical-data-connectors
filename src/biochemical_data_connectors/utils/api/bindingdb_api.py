import time
import requests
import logging
from typing import List, Optional, Any

from biochemical_data_connectors.utils.api.base_api import BaseApiClient
from biochemical_data_connectors.constants import RestApiEndpoints


class BindingDbApiClient(BaseApiClient):
    """
    A client for interacting with the BindingDB REST API.

    This class encapsulates direct communication with BindingDB endpoints,
    using the persistent session and retry logic inherited from BaseAPIClient.
    """
    def __init__(self, logger: Optional[logging.Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_actives_from_target_uniprot(
        self,
        uniprot_id: str,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,  # In nM.
    ) -> List[Any]:
        """
        Queries BindingDB for all affinity data for a given UniProt ID.

        This method constructs the appropriate URL for the `getLigandsByUniprot`
        endpoint, applies the affinity cutoff at the API level, and then performs
        a local filter for the desired bioactivity measure types.

        Parameters
        ----------
        uniprot_id : str
            The UniProt accession ID of the target.
        bioactivity_measures : List[str]
            A list of activity types to keep (e.g., ['Ki', 'IC50']).
        bioactivity_threshold : float, optional
            The affinity cutoff in nM to be passed to the API. If None, no
            cutoff is applied.

        Returns
        -------
        List[Dict]
            A list of dictionaries, where each dictionary is an affinity record
            from the BindingDB API that matches the filter criteria.
        """
        bindingdb_start = time.time()

        cutoff_str = f';{int(bioactivity_threshold)}' if bioactivity_threshold is not None else ''
        url = RestApiEndpoints.BINDINGDB_LIGANDS_FROM_UNIPROT_ID.url(
            uniprot_id=uniprot_id,
            cutoff_str=cutoff_str
        )
        self._logger.info(f'Querying BindingDB for target: {uniprot_id}')
        try:
            response = self._session.get(url)
            response.raise_for_status()
            data = response.json()
            bdb_affinities = data.get('getLindsByUniprotResponse', {}).get('bdb.affinities', [])

            bindingdb_end = time.time()
            self._logger.info(f'BindingDB total query time: {round(bindingdb_end - bindingdb_start)} seconds')
            if not bdb_affinities:
                self._logger.warning(f'No BindingDB actives found for {uniprot_id}.')
                return []

            filtered_bdb_actives = [
                record for record in bdb_affinities if record.get('bdb.affinity_type') in bioactivity_measures
            ]

            return filtered_bdb_actives

        except requests.exceptions.RequestException as e:
            self._logger.error(f'Error querying BindingDB: {e}')
            return []
