import time
from logging import Logger

import requests
import logging
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from typing import List, Dict, Optional

from chembl_webresource_client.new_client import new_client

from biochemical_data_connectors.utils.api.base_api import BaseAPIClient
from biochemical_data_connectors.constants import RestApiEndpoints

class ChEMBLAPIClient(BaseAPIClient):
    def __init__(self, logger: Optional[Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_activities_for_target(self, target_chembl_id: str, activity_types: List[str]) -> List[Dict]:
        """Paginates through the ChEMBL activity API to fetch all records for a target."""
        all_records = []
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type__in": ",".join(activity_types),  # Use '__in' for multiple values
        }

        limit = 1000
        offset = 0
        chembl_activity_url = RestApiEndpoints.CHEMBL_ACTIVITY.url()

        while True:
            page_params = {
                **params,
                "limit": limit,
                "offset": offset
            }
            try:
                response = self._session.get(chembl_activity_url, params=page_params, timeout=15)
                response.raise_for_status()
                data = response.json()
                records = data.get('activities', [])
                if not records:
                    break
                all_records.extend(records)
                offset += limit
            except requests.exceptions.RequestException as e:
                self._logger.error(f"ChEMBL API request failed for target {target_chembl_id}: {e}")
                break

        chembl_end = time.time()
        self._logger.info(f'ChEMBL total query time: {round(chembl_end - chembl_start)} seconds')

        return all_records
