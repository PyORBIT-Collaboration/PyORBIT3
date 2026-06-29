def multiDimDoubleArray(*dims: int) -> list[float]:
    """
    Creates multi-dimensional arrays with doubles, such as a[i][k][j].

    Parameters
    ----------
    dims : int

    Returns
    -------
    list[float]

    Note
    ----
    All elements are initialized to 0.

    Examples
    --------
    >>> a = multiDimArray(5,10,2)
    >>> a = multiDimArray(*[5,10,2])  # equivalent
    >>> a[1][2][1]
    0.0
    """
    res = []
    if len(dims) == 1:
        for j in range(dims[0]):
            res.append(0.0)
    else:
        dims_rest = dims[1 : len(dims)]
        for j in range(dims[0]):
            res.append(multiDimDoubleArray(*dims_rest))
    return res


def multiDimIntArray(*dims: int) -> list[int]:
    """
    Creates multi-dimensional arrays with integers, such as a[i][k][j].

    Parameters
    ----------
    dims : int

    Returns
    -------
    list[int]

    Examples
    --------
    >>> a = multiDimArray(5,10,2)
    >>> a = multiDimArray(*[5,10,2])  # equivalent
    >>> a[1][2][1]
    0
    """
    res = []
    if len(dims) == 1:
        for j in range(dims[0]):
            res.append(0)
    else:
        dims_rest = dims[1 : len(dims)]
        for j in range(dims[0]):
            res.append(multiDimIntArray(*dims_rest))
    return res
