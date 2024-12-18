# quangpy: Python package for Quantitative Genetics Analysis

## Overview

`quangpy` mainly provides Joint Scale Analysis and Variance Decomposition Analysis.

## Installation

```bash
pip install quangpy

```

## Getting Started

### input data format

prepare data should be sep by tab, as the following format with first column
The input file should sep by tab and provide the following columns with the same order:

- population: population name
- mean: mean of the population phenotype
- var: variance

The following is an example of the input data format:

|     | mean  | var                  |
| --- | ----- | -------------------- |
| b1l | 7.5   | 0.010000000000000002 |
| b1s | 5.575 | 0.004096             |
| f1  | 5.5   | 0.007395999999999999 |
| f2  | 6.595 | 0.013923999999999999 |
| p1  | 9.125 | 0.008281             |
| p2  | 5.475 | 0.003249             |

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
