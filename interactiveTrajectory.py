!pip install scanpy
!pip install plotly
import scanpy as sc
import plotly.graph_objects as go

## Load the data
adata = sc.datasets.paul15()

# preprocessing

sc.pp.recipe_zheng17(adata)

# 2D
def plot_interactive_trajectory_2d(adata):

    # Compute neighborhood graph
    sc.pp.neighbors(adata)
    
    # Perform diffusion map embedding
    sc.tl.diffmap(adata)
    
    # Perform trajectory inference
    sc.tl.dpt(adata)
    
    # Create an interactive scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=adata.obsm['X_diffmap'][:, 0], y=adata.obsm['X_diffmap'][:, 1], mode='markers', marker=dict(color=adata.obs['dpt_pseudotime'],colorscale='viridis'),text=adata.obs_names))
    
    # Add layout and axis labels
    fig.update_layout(title='Interactive Trajectory Plot', xaxis_title='Diffmap Component 1',yaxis_title='Diffmap Component 2')
    
    # Show the plot
    fig.show(config={"displayModeBar": True}, auto_open=True)


# 3D
def plot_interactive_trajectory_3d(adata):
    #Plot interactive 3D trajectory visualization.
    # Compute neighborhood graph
    sc.pp.neighbors(adata)
    
    # Perform diffusion map embedding
    sc.tl.diffmap(adata)
    
    # Perform trajectory inference (e.g., pseudotime estimation)
    sc.tl.dpt(adata)
    
    # Create an interactive 3D scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=adata.obsm['X_diffmap'][:, 0], y=adata.obsm['X_diffmap'][:, 1], z=adata.obsm['X_diffmap'][:, 2],mode='markers',marker=dict(color=adata.obs['dpt_pseudotime'],colorscale='viridis'),text=adata.obs_names))
    
    # Add layout and axis labels
    fig.update_layout(title='Interactive 3D Trajectory Plot',scene=dict(xaxis_title='Diffmap Component 1', yaxis_title='Diffmap Component 2',zaxis_title='Diffmap Component 3'))
    
    # Show the plot
    fig.show(config={"displayModeBar": True}, auto_open=True)

# testing of the 2 previous functions    
plot_interactive_trajectory_2d(adata)    
plot_interactive_trajectory_3d(adata)   
