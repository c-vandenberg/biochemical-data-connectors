import logging
from abc import ABC, abstractmethod
from typing import List, Optional

from biochemical_data_connectors.models import BioactiveCompound


class BaseBioactivesConnector(ABC):
    """
    Abstract base class for extracting bioactive compounds from a data source.

    Attributes
    ----------
    _bioactivity_measures : List[str]
        A prioritized list of bioactivity measurement types to filter on
        (e.g., ['Kd', 'Ki', 'IC50']).
    _bioactivity_threshold : float, optional
        The maximum potency value (in nM) to consider a compound bioactive.
    _cache_dir : str
        The local directory path for storing cached data.
    _logger : logging.Logger
        A logger instance for logging messages.
    """
    def __init__(
        self,
        bioactivity_measures: List[str],
        bioactivity_threshold: Optional[float] = None,
        cache_dir: str = './data/cache',
        logger: Optional[logging.Logger] = None
    ):
        self._bioactivity_measures = bioactivity_measures
        self._bioactivity_threshold = bioactivity_threshold
        self._cache_dir = cache_dir
        self._logger = logger if logger else logging.getLogger(__name__)

    @abstractmethod
    def get_bioactive_compounds(self, target_uniprot_id: str, force_refresh: bool = False) -> List[BioactiveCompound]:
        """
        Retrieve a list of canonical SMILES for bioactive compounds for a given target.

        Parameters
        ----------
        target_uniprot_id : str
            The target identifier (UniProt accession, e.g. "P00533").
        force_refresh : bool, optional
            If True, ignores any existing cache and forces a new API call.
            Default is False.

        Returns
        -------
        List[BioactiveCompound]
            A list of structured BioactiveCompound objects.
        """
        pass
