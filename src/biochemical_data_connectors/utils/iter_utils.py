from typing import Iterator, List, Sequence


def batch_iterable(iterable: List, n: int = 1000) -> Iterator[List]:
    """
    Splits a sequence into smaller, successive batches of a given size.

    Parameters
    ----------
    iterable : Sequence[Any]
        The sequence (e.g., list) to split into batches.
    n : int, optional
        The desired batch size. Default is 1000.

    Yields
    ------
    Iterator[Sequence[Any]]
        A generator that yields one batch (as a list or slice) at a time.

    Examples
    --------
    >>> my_list = [1, 2, 3, 4, 5, 6, 7]
    >>> for batch in batch_iterable(my_list, 3):
    ...     print(list(batch))
    [1, 2, 3]
    [4, 5, 6]
    [7]
    """
    iter_len = len(iterable)
    for idx in range(0, iter_len, n):
        yield iterable[idx:idx + n]
