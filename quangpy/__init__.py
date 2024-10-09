"""Simple package for quantitative genetic analysis."""

from .joint_scale import JointScale
from .var_decomp import VarianceDecomposition
from .wls import WeightedLeastSquare

__all__ = ["JointScale", "VarianceDecomposition", "WeightedLeastSquare"]
