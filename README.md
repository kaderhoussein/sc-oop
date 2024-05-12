
# Interactive diffusion map for Single-Cell Data Analysis

This repository contains adapted functions from the scanpy package and custom plots to visualize single-cell data. The repository includes a requirements.txt file with the necessary libraries to run the code.

## Getting Started

### Prerequisites

Before starting, make sure to clone this repository

```bash
git clone https://github.com/kaderhoussein/sc_oop.git
```
You will also need to install the required libraries either using pip or directly using the requirements.txt file

```bash
cd sc_oop
pip install -r requirements.txt
```

## Importing the Package
To use the interactive_diffmap_2D function, you need to import the following modules:

(To use this option, make sure you are placed in this repository root folder.)

```python
import scanpy as sc
import plotly.graph_objects as go
from scoop.plotting.interactive_plots.iscatterplots import interactive_diffmap_2D, interactive_diffmap_3D, interactive_UMAP_2D
```


## Example of use

Here's an example of how to use the interactive_diffmap_2D and the interactive_UMAP_2D  functions:

```python
adata = sc.datasets.paul15()
interactive_diffmap_2D(adata)
interactive_UMAP_2D(adata)
```
## Install

```bash
# From sc_oop root
pip install -e .
```

Then, you could run the following:

```python
import scanpy as sc
import plotly.graph_objects as go
import scoop.plotting.interactive_plots.iscatterplots

adata = sc.datasets.paul15()
scoop.plotting.interactive_plots.iscatterplots.interactive_diffmap_2D(adata)
```

## Authors
 - Benamad Kader Houssein
 - Abdou Khadre Djeylani Diouf
