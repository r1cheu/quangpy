"""Joint Scale Analysis."""

from pathlib import Path

import numpy as np
from scipy.stats import chi2

from .utils import print_list
from .wls import WeightedLeastSquare


class JointScale:
    """
    Joint Scale Analysis.

    Parameters
    ----------
    input_data : str | Path
        Input data file path.

    See Also
    --------
    quantpy.WeightedLeastSquare : Weighted Least Square.
    quantpy.VarianceDecomposition : Variance Decomposition Analysis.

    Notes
    -----
    The input file should sep by tab and provide the following columns with the same order:

    - population: population name
    - mean: mean of the population phenotype
    - var: variance

    Examples
    --------
    >>> from quantpy import JointScale
    >>> input_data = "path/to/input/data.txt"
    >>> effect = ("m", "a",)
    >>> joint_scale = JointScale(input_data)
    >>> joint_scale.fit(effect)
    =================================
          coef  std err   P>|t|
    =================================
       m    6.6740   0.3187  0.0000
       a    1.6333   0.4398  0.0206
    chi2  345.0712           0.0000
    =================================
    """

    def __init__(self, input_data: str | Path) -> None:
        self._wls = WeightedLeastSquare(input_data)

    def _check_effect(self, effect: tuple[str]) -> tuple[str]:
        """
        Check the effect is valid.

        Only "m", "a", "d", "aa", "ad", "dd" are allowed. And the number of effect
        should be less than the df of the model(5 here).

        Parameters
        ----------
        effect : tuple[str]
            The effect to check.

        Returns
        -------
        tuple[str]
            Valid effect.
        """
        all_effect = ["m", "a", "d", "aa", "ad", "dd"]
        effect = list(effect)
        if "m" not in effect:
            effect.insert(0, "m")
        if not set(effect).issubset(set(all_effect)):
            msg = f"Invalid effect, only {all_effect} are allowed."
            raise ValueError(msg)
        if len(effect) > 5:
            msg = (
                "The number of effect should be less than the df of the model(5 here)."
                f" but got {len(effect)} effect"
            )
            raise ValueError(msg)
        return effect

    def fit(self, effect: tuple[str]) -> None:
        """
        Fit the model with the given effect.

        Parameters
        ----------
        effect : tuple[str]
            Effect to fit the model.
        """
        effect = self._check_effect(effect)
        result = self._wls.fit(effect)
        self._format_result(result)

    def _format_result(self, result) -> None:
        """
        Format the result to readable table.

        Parameters
        ----------
        result : statsmodels.regression.linear_model.RegressionResultsWrapper
            Result of the statsmodels regression.
        """
        lines = [["", "coef", "std err", "P>|t|"]]  # header
        for i in result.model.exog_names:  # for each effect
            sig = (
                "**"
                if result.pvalues[i] < 0.01
                else "*"
                if result.pvalues[i] < 0.05
                else ""
            )
            lines.append(
                [
                    i,
                    f"{result.params[i]:.4f}",
                    f"{result.bse[i]:.4f}",
                    f"{result.pvalues[i]:.4f}{sig}",
                ]
            )
        chi_result = [f"{i:.4f}" for i in self._chi_prob(result)]
        chi_result.insert(1, "")
        lines.append(["chi2", *chi_result])
        print_list(lines)

    def _chi_prob(self, result) -> tuple[float, float]:
        """
        Calculate chi_squared value and the probability of chi_squared.

        Parameters
        ----------
        result : statsmodels.regression.linear_model.RegressionResultsWrapper
            Result of the statsmodels regression.

        Returns
        -------
        tuple[float, float]
            Chi_squared value and the probability of chi_squared.
        """
        residuals = result.resid
        chi_squared = np.sum((residuals**2) * result.model.weights)
        return chi_squared, chi2.sf(chi_squared, result.df_resid)
