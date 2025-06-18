import os
import json
import pickle
import logging
from datetime import datetime
from typing import Optional, Any, Callable


def get_cached_or_fetch(
    cache_file_path: str,
    fetch_function: Callable[[], Any],
    data_type: str = 'bioactive',
    use_pickle: bool = False,
    force_refresh: bool = False,
    logger: Optional[logging.Logger] = None
) -> Optional[Any]:
    """
    Generic method to retrieve data from a cache file or by executing a
    data fetching function if the cache does not exist or a refresh is forced.

    Parameters
    ----------
    cache_file_path : str
        The path to the cache file.
    fetch_function : Callable[[], Any]
        A no-argument function that will be called to get fresh data if needed.
    data_type : str, optional
        A descriptive name for the data being processed, used for logging.
        Default is 'items'.
    use_pickle : bool, optional
        If True, use the binary pickle format for serialization. If False
        (default), use JSON.
    force_refresh : bool, optional
        If True, ignores the cache and always executes the fetch_function.
    logger : logging.Logger, optional
        A logger instance for status messages. If `None` (default), messages
        are printed to standard output.

    Returns
    -------
    Any
        The data from the cache or the fetch_function.
    """
    cache_is_valid = False
    data = None
    # 1) If manual `--force-refresh` flag is True, and cache file exists, load data from
    #    cache file.
    if not force_refresh and os.path.exists(cache_file_path):
        if use_pickle:
            with open(cache_file_path, 'rb') as f:  # 'rb' for read binary
                cache_data = pickle.load(f)
        else:
            with open(cache_file_path, 'r') as f:
                cache_data = json.load(f)

        cache_message = f"Found valid cache at {cache_file_path}. Loading {data_type} data from file."
        logger.info(cache_message) if logger else print(cache_message)

        data = cache_data["data"]
        cache_is_valid = True

    # 2) If cache file doesn't exist or `force_refresh == True`, call the provided fetch_function.
    if not cache_is_valid:
        pre_api_message = f"Fetching fresh {data_type} data from API..."
        logger.info(pre_api_message) if logger else print(pre_api_message)

        data = fetch_function()

        # 3) Save the new results to the cache file with a current timestamp.
        if data:
            cache_content = {
                "timestamp": datetime.now().isoformat(),
                "data": data
            }
            os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
            if use_pickle:
                with open(cache_file_path, 'wb') as f:
                    pickle.dump(cache_content, f)
            else:
                with open(cache_file_path, 'w') as f:
                    json.dump(cache_content, f, indent=4)

            post_api_message = f"Saved {len(data)} {data_type} items to cache file: {cache_file_path}"
            logger.info(post_api_message) if logger else print(post_api_message)

    return data