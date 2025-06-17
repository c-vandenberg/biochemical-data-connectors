import time
from logging import Logger

import requests
import logging
from typing import List, Dict, Optional

from biochemical_data_connectors.utils.api.base_api import BaseAPIClient
from biochemical_data_connectors.constants import RestApiEndpoints


class ChEMBLAPIClient(BaseAPIClient):
    """
    A client for interacting with the ChEMBL API.

    This class encapsulates all direct interactions with the ChEMBL REST API.
    A persistent `requests.Session` instance with a retry strategy handles
    transient network errors and API rate-limiting.
    """
    def __init__(self, logger: Optional[Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_activities_for_target(self, target_chembl_id: str, activity_types: List[str]) -> List[Dict]:
        """
        Paginates through the ChEMBL activity API to fetch all records.

        Parameters
        ----------
        target_chembl_id : str
            The ChEMBL ID of the target to query for activities.
        activity_types : List[str]
            A list of standard activity types to filter for (e.g., ['Kd', 'Ki']).

        Returns
        -------
        List[Dict]
            A list of dictionary objects, where each dictionary is a full
            activity record from the ChEMBL API.
        """
        # 1) Initialize a list to store all fetched records and start the timer.
        chembl_start = time.time()
        all_records = []

        # 2) Build the base parameters for ChEMBL REST API.
        #    The '__in' suffix allows filtering by multiple standard_type values.
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type__in": ",".join(activity_types),
        }

        # 3) Set up variables for paginating through the API results.
        limit = 1000 # Number of records to fetch per page
        offset = 0 # Starting point for the record set
        chembl_activity_url = RestApiEndpoints.CHEMBL_ACTIVITY.url()

        # 4) Loop continuously to fetch all pages of data until no more records are returned.
        while True:
            # 4.1) Set the parameters for the current page, including limit and offset.
            page_params = {
                **params,
                "limit": limit,
                "offset": offset
            }
            try:
                # 4.2) Make the GET request using the persistent session.
                response = self._session.get(chembl_activity_url, params=page_params, timeout=15)
                response.raise_for_status()
                data = response.json()

                # 4.3) Parse the JSON response and extract the list of activities.
                records = data.get('activities', [])

                # 4.4) If no records are returned, it's the last page, so break the loop.
                #      Otherwise, add the fetched records and increment the offset to get
                #      the next page in the iteration.
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
