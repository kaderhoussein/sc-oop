# Single-Cell RNA Sequencing Trajectory Visualization
This project provides Python scripts for visualizing trajectory inferences in single-cell RNA sequencing (scRNA-seq) data using diffusion maps and trajectory inference techniques. The visualization is interactive and can be viewed in both 2D and 3D.

# Usage
Clone or download this repository to your local machine.
Run the Python script trajectory_visualization.py with your scRNA-seq data. Make sure to replace adata with your own dataset.
The script will preprocess the data, perform trajectory inference, and generate interactive plots in both 2D and 3D.

python trajectory_visualization.py

# Example
Here's an example of how to use the provided functions:

## Import necessary libraries

`import scanpy as sc`

`import plotly as pl`

`import plotly.graph_objects as go`

## Load example scRNA-seq data

`adata = sc.datasets.paul15()`

## Run 2D trajectory visualization

`plot_interactive_trajectory_2d(adata)`

## Run 3D trajectory visualization

`plot_interactive_trajectory_3d(adata)`

## Team members
- Benamad Kader Houssein
- Abdou Khadre DIOUF
