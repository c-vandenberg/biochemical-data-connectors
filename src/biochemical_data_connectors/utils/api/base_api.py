import requests
import logging
from requests.adapters import HTTPAdapter, Retry
from typing import Optional


class BaseAPIClient:
    """
    A base client for API interactions with a persistent, robust session.

    Attributes
    ----------
    _logger : logging.Logger
        A logger instance for logging messages.
    _session : requests.Session
        A configured requests session with automatic retries for server errors.
    """
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
            status_forcelist=[429, 500, 502, 503, 504]
        )
        # Mount the retry strategy to the session for all HTTPS requests.
        session.mount("https://", HTTPAdapter(max_retries=retries))

        return session