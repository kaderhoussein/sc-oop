# Single-Cell RNA Sequencing Trajectory Visualization
This project provides Python scripts for visualizing trajectory inferences in single-cell RNA sequencing (scRNA-seq) data using diffusion maps and trajectory inference techniques. The visualization is interactive and can be viewed in both 2D and 3D.

# Installation
To run the scripts, you need to have Python installed on your system along with the following packages:

scanpy: For preprocessing and analyzing scRNA-seq data.
plotly: For creating interactive plots.
You can install these dependencies using pip:

`pip install scanpy`

`pip install plotly`

# Usage
Clone or download this repository to your local machine.
Run the Python script trajectory_visualization.py with your scRNA-seq data. Make sure to replace adata with your own dataset.
The script will preprocess the data, perform trajectory inference, and generate interactive plots in both 2D and 3D.

python trajectory_visualization.py

# Example
Here's an example of how to use the provided functions:

## Import necessary libraries

`import scanpy as sc`

## Load example scRNA-seq data

`adata = sc.datasets.paul15()`

## Run 2D trajectory visualization

`plot_interactive_trajectory_2d(adata)`

## Run 3D trajectory visualization

`plot_interactive_trajectory_3d(adata)`
