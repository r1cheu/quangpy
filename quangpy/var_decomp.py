"""Variance decomposition analysis."""

from functools import cached_property

import numpy as np
from scipy.stats import f

from .utils import print_list, tofloat
from .wls import WeightedLeastSquare


class VarianceDecomposition:
    """
    Perform variance decomposition analysis.

    Parameters
    ----------
    input_data : str
        The path to the input CSV file containing the data.

    See Also
    --------
    quantpy.WeightedLeastSquare : Weighted Least Square.
    quantpy.JointScale : Joint Scale Analysis.

    Notes
    -----
    The input file should sep by tab and provide the following columns with the same order:

    - population: population name
    - mean: mean of the population phenotype
    - var: variance

    Examples
    --------
    >>> from quantpy import VarianceDecomposition
    >>> input_data = "path/to/input/data.txt"
    >>> effect = ("m", "a", "d", "aa", "dd")
    >>> vardecomp = VarianceDecomposition(input_data)
    >>> vardecomp.fit(effect)
    ===========================================================================================
               MS  Prob(F)       MS  Prob(F)       MS  Prob(F)       MS  Prob(F)    propotion
    ===========================================================================================
     [a]  1189.96   0.0206  1189.96   0.0001  1189.96   0.0003  1189.96   0.0142       0.7755
     [d]                     340.70   0.0006   340.70   0.0010   340.70   0.0265       0.2220
    [aa]                                         3.72   0.0774     3.72   0.2411       0.0024
    [dd]                                                           0.06   0.8012       0.0000
      R2   0.7752            0.9972            0.9996            0.9996           SSy=1535.04
    ===========================================================================================
    """

    def __init__(self, input_data: str) -> None:
        self._vardecomp = VarDecomp(input_data)

    def _check_effect(self, effect):
        """
        Check if the provided effects are valid.

        Only "a", "d", "aa", "dd", "ad" are allowed. The number of effects should be equal to 4.

        Parameters
        ----------
        effect : list
            A list of effects to be checked.

        Returns
        -------
        list
            The validated list of effects.

        Raises
        ------
        ValueError
            If the effects are not valid or the number of effects is not equal to 4.
        """
        allow_eff = ["a", "d", "aa", "ad", "dd"]
        effect = list(effect)
        if not set(effect).issubset(set(allow_eff)):
            msg = "Invalid effect, only `a`, `d`, `aa`, `ad`, `dd` are allowed."
            raise ValueError(msg)
        if len(effect) > 5:
            msg = (
                "The number of effect should be less than the df of the model(5 here)."
                f" but got {len(effect)} effect"
            )
            raise ValueError(msg)
        return effect

    def fit(self, effect):
        """
        Fit the variance decomposition model using the specified effects.

        Parameters
        ----------
        effect : list
            A list of effects to be used in the model.
        """
        effect = self._check_effect(effect)
        _effect = ["m"]
        SS_effect = []
        f_prob = []
        for idx, e in enumerate(effect):
            prev_SSR = self._vardecomp.SSR if idx > 0 else 0
            _effect.append(e)
            _ = self._vardecomp.fit(_effect)
            SS_effect.append(self._vardecomp.SSR - prev_SSR)
            f_prob.append(self._vardecomp.f_prob(SS_effect))

        self._format_result(effect, SS_effect, f_prob)

    def _format_result(self, effect, SS_effect, f_prob):
        """
        Format result to readable table.

        Parameters
        ----------
        effect : list
            A list of effects used in the model.
        SS_effect : list
            A list of sum of squares for each effect.
        f_prob : list
            A list of F-test probabilities for each effect.
        """
        n = len(effect)
        lines = [[""] + ["MS", "Prob(F)"] * n + ["propotion"]]
        propotion = [SS_effect[i] / self._vardecomp.SSy for i in range(n)]
        r2 = [sum(SS_effect[: i + 1]) / self._vardecomp.SSy for i in range(n)]

        # handle effect line
        for idx, eff in enumerate(effect):
            _line = [f"[{eff}]"]
            for _ in range(idx):
                _line.extend([""] * 2)  # add empty cells

            for p in f_prob[idx:]:
                sig = "**" if p[idx] < 0.01 else "*" if p[idx] < 0.05 else ""
                _line.extend(
                    [f"{SS_effect[idx]:.2f}", f"{p[idx]:.4f}{sig}"]
                )  # add SS and p-value

            _line.append(f"{propotion[idx]:.4f}")  # add propotion
            lines.append(_line)

        # handle r2 line
        r2_line = ["R2"]
        for _r2 in r2:
            r2_line.extend([f"{_r2:.4f}", ""])
        r2_line.append(f"SSy={self._vardecomp.SSy:.2f}")
        lines.append(r2_line)

        print_list(lines)


class VarDecomp(WeightedLeastSquare):
    """
    Weighted Least Square for variance decomposition.

    Provide some additional methods to calculate the sum of squares and F-test probabilities.
    This class is not intended to be used directly.

    Parameters
    ----------
    input_data : str
        The path to the input file containing the data.
    """

    def __init__(self, input_data: str) -> None:
        super().__init__(input_data)

    @cached_property
    @tofloat
    def SSy(self):
        """
        Calculate the sum of squares of the response variable.

        Returns
        -------
        float
            The sum of squares of the response variable.
        """
        return self._y.T @ self._lambda @ self._y - self._correction

    @cached_property
    @tofloat
    def _correction(self):  # numpydoc ignore=GL08
        return np.sum(self._weight * self._y) ** 2 / np.sum(self._weight)

    @property
    @tofloat
    def SSR(self):
        """
        Calculate the regression sum of squares.

        Returns
        -------
        float
            The regression sum of squares.

        Raises
        ------
        ValueError
            If the model has not been fitted yet.
        """
        if self.result is None:
            msg = "Please fit the model first"
            raise ValueError(msg)
        result = self.result
        return (
            result.params @ result.model.exog.T @ self._lambda @ result.model.endog
            - self._correction
        )

    @property
    @tofloat
    def _SSr(self):  # numpydoc ignore=GL08
        if self.result is None:
            msg = "Please fit the model first"
            raise ValueError(msg)
        return self.SSy - self.SSR

    @property
    @tofloat
    def _df(self):  # numpydoc ignore=GL08
        if self.result is None:
            msg = "Please fit the model first."
            raise ValueError(msg)
        return self.result.df_resid

    @property
    @tofloat
    def _MSr(self):  # numpydoc ignore=GL08
        return self._SSr / self._df

    def f_prob(self, SSR):
        """
        Calculate the F-test probabilities for the given SSR values.

        Parameters
        ----------
        SSR : list
            A list of SSR values for which to calculate the F-test probabilities.

        Returns
        -------
        list
            A list of F-test probabilities corresponding to the given SSR values.
        """

        @tofloat
        def f_prob(ssr):  # numpydoc ignore=GL08
            return f.sf(ssr / self._MSr, 1, self._df)

        return [f_prob(ssr) for ssr in SSR]
