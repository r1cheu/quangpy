"""Weighted Least Squares (WLS) regression module."""

import numpy as np
import pandas as pd
import statsmodels.api as sm


class WeightedLeastSquare:
    """
    Perform Weighted Least Squares (WLS) regression.

    Parameters
    ----------
    input_data : str
        The path to the input file.

    Attributes
    ----------
    DM : pd.DataFrame
        A predefined design matrix used for the regression.
    """

    DM = pd.DataFrame(
        [
            [1.0, 1.0, 0.0, 1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.5, 0.0, 0.0, 0.25],
            [1.0, 0.5, 0.5, 0.25, 0.25, 0.25],
            [1.0, -0.5, 0.5, 0.25, -0.25, 0.25],
        ],
        index=pd.Index(["p1", "p2", "f1", "f2", "b1l", "b1s"]),
        columns=pd.Index(["m", "a", "d", "aa", "ad", "dd"]),
    )

    def __init__(self, input_data: str) -> None:
        """
        Initialize the WeightedLeastSquare with input data from a CSV file.

        Parameters
        ----------
        input_data : str
            The path to the input CSV file containing the data.
        """
        self._data = pd.read_csv(input_data, index_col=0, sep="\t")
        self._weight = 1 / self._data["var"].to_numpy()
        self._lambda = np.diag(self._weight)

        self._y = self._data["y"].to_numpy()

    def fit(self, effect: tuple[str]):
        """
        Fit the WLS model using the specified effect and return the result.

        Parameters
        ----------
        effect : str
            The effect to be used in the regression model.

        Returns
        -------
        sm.regression.linear_model.RegressionResultsWrapper
            The result of the fitted WLS model.
        """
        X = self.DM.loc[self._data.index, effect]
        model = sm.WLS(self._y, X, weights=self._weight)
        result = model.fit()
        self.result = result
        return result
