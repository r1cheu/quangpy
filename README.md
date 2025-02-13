# quangpy: a python package for **q**uantitative **g**enetics analysis

[![](https://img.shields.io/pypi/v/quangpy.svg?color=brightgreen)](https://pypi.org/pypi/quangpy/)
[![](https://img.shields.io/github/issues/r1cheu/quangpy?color=green)](https://github.com/r1cheu/quangpy/issues/new)
[![](https://img.shields.io/github/license/r1cheu/quangpy)]()

`quangpy` mainly provides Joint Scale Analysis and Variance Decomposition Analysis.

## Installation

```bash
pip install quangpy

```

## Getting Started

### input data format

The input file should be separated by tab and provide the following columns in the same order:

- population: population name
- mean: mean of the population phenotype
- var: variance of the population phenotype

The following is an example of the input data format:

|     | mean  | var      |
| --- | ----- | -------- |
| b1l | 7.5   | 0.010000 |
| b1s | 5.575 | 0.004096 |
| f1  | 5.5   | 0.007396 |
| f2  | 6.595 | 0.013924 |
| p1  | 9.125 | 0.008281 |
| p2  | 5.475 | 0.003249 |

### Joint Scale Analysis

```python
from quangpy import JointScale

input_data = "/path/to/input/data"

js = JointScale(input_data)
js.fit(("m", "a"))

```

### Variance Decomposition Analysis

```python
from quangpy import VarianceDecomposition

input_data = "/path/to/input/data"
vd = VarianceDecomposition(input_data)
vd.fit(("m", "a", "d", "aa", "dd"))
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
