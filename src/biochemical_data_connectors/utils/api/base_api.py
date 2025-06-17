import requests
import logging
from requests.adapters import HTTPAdapter, Retry
from typing import List, Dict, Optional


class BaseAPIClient:
    def __init__(self, logger: Optional[logging.Logger] = None):
        self._logger = logger if logger else logging.getLogger(__name__)
        self._session = self._create_session()

    @staticmethod
    def _create_session() -> requests.Session:
        """
        Create a `requests.Session` instance with a robust retry strategy.

        The retry strategy automatically handles transient network issues and
        common server-side errors (like 5xx status codes), making the client
        more resilient.

        Returns
        -------
        requests.Session
            A configured session object with mounted retry logic.
        """
        session = requests.Session()
        retries = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[500, 502, 503, 504]
        )
        # Mount the retry strategy to the session for all HTTPS requests.
        session.mount("https://", HTTPAdapter(max_retries=retries))

        return session