"""Provide utility functions for the project."""


def tofloat(func):
    """
    Decorator that converts the result of a function to a float.

    Parameters
    ----------
    func : callable
        The function to be decorated.

    Returns
    -------
    callable
        The decorated function that returns a float.
    """

    def wrapper(*args, **kwargs):  # numpydoc ignore=GL08
        return float(func(*args, **kwargs))

    return wrapper


def print_list(data):
    """
    Print a 2D list in a formatted table with columns aligned.

    Parameters
    ----------
        data : list[list[str]]
            The 2D list to be printed.
    """

    def calc_max_width(data):  # numpydoc ignore=GL08
        max_col_width = []
        for row in data:
            for idx, col in enumerate(row):
                if len(max_col_width) <= idx:
                    max_col_width.append(0)
                max_col_width[idx] = max(max_col_width[idx], len(str(col)))
        return max_col_width

    def print_split_line(max_col_width, element="="):  # numpydoc ignore=GL08
        print("".join([element * (col + 2) for col in max_col_width]))  # noqa:T201

    max_col_width = calc_max_width(data)
    print_split_line(max_col_width, "=")

    for out_idx, row in enumerate(data):
        for idx, col in enumerate(row):
            print(f"{col:<{max_col_width[idx]}} ", end=" ")  # noqa:T201
        print()  # noqa:T201
        if out_idx == 0:
            print_split_line(max_col_width, "=")

    print_split_line(max_col_width, "=")
