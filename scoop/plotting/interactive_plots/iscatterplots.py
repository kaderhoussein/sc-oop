import scanpy as sc
import plotly.graph_objects as go

# 2D

def interactive_diffmap_2D(adata):

    # preprocessing using the log1p method
    sc.pp.log1p(adata)   

    # Compute neighborhood graph
    sc.pp.neighbors(adata)

    # Perform diffusion map embedding
    sc.tl.diffmap(adata)

    # Perform trajectory inference
    sc.tl.dpt(adata)

    # Create the interactive 2D scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=adata.obsm['X_diffmap'][:, 0],
                             y=adata.obsm['X_diffmap'][:, 1],
                             mode='markers', 
                             marker=dict(color=adata.obs.iloc[:,1], colorscale='viridis', colorbar=dict(title='Pseudotime')), #2nd column obs
                             text=adata.obs.iloc[:,0],#1st column obs
                             #customdata=adata.obs['n_counts_all'],
                             hovertemplate='Pseudotime: %{marker.color:.2f}<br>Cluster: %{text}'
                             )
                 )

    # Add layout and axis labels
    fig.update_layout(title='Interactive Trajectory Plot',
                      xaxis_title='Diffmap Component 1',
                      yaxis_title='Diffmap Component 2')

    # Add buttons for selecting method
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                buttons=list([
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:, 0]],
                               "y": [adata.obsm['X_diffmap'][:, 1]],
                               "marker.color": [adata.obs.iloc[:,1]],
                               "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 1",
                        method="restyle"
                    ),
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:,1]],
                               "y": [adata.obsm['X_diffmap'][:, 2]],
                               "marker.color": [adata.obs.iloc[:,1]],
                                "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 2",    
                        method="restyle"
                    ),
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:,2]],
                               "y": [adata.obsm['X_diffmap'][:, 3]],
                               "marker.color": [adata.obs.iloc[:,1]],
                                "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 3",    
                        method="restyle"
                    
                    )
                ])
            )
        ]
    )

    # Show the plot
    plot = fig.show(config={"displayModeBar": True}, auto_open=True)
    return(plot)


# 3D

def interactive_diffmap_3D(adata):

    # preprocessing using the log1p method
    sc.pp.log1p(adata)   

    # Compute neighborhood graph
    sc.pp.neighbors(adata)

    # Perform diffusion map embedding
    sc.tl.diffmap(adata)

    # Perform trajectory inference
    sc.tl.dpt(adata)

    # Create the interactive 3D scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=adata.obsm['X_diffmap'][:, 0],
                             y=adata.obsm['X_diffmap'][:, 1],
                             z=adata.obsm['X_diffmap'][:, 2],
                             mode='markers', 
                             marker=dict(color=adata.obs.iloc[:,1], colorscale='viridis', colorbar=dict(title='Pseudotime')), #2nd column obs
                             text=adata.obs.iloc[:,0],#1st column obs
                             #customdata=adata.obs['n_counts_all'],
                             hovertemplate='Pseudotime: %{marker.color:.2f}<br>Cluster: %{text}'
                             )
                 )

    # Add layout and axis labels
    fig.update_layout(title='Interactive Trajectory Plot',
                      scene=dict(xaxis_title='Diffmap Component 1',
                                 yaxis_title='Diffmap Component 2',
                                 zaxis_title='Diffmap Component 3'))

    # Add buttons for selecting method
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                buttons=list([
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:, 0]],
                               "y": [adata.obsm['X_diffmap'][:, 1]],
                               "z": [adata.obsm['X_diffmap'][:, 2]],
                               "marker.color": [adata.obs.iloc[:,1]],
                               "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 1",
                        method="restyle"
                    ),
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:,1]],
                               "y": [adata.obsm['X_diffmap'][:, 2]],
                               "z": [adata.obsm['X_diffmap'][:, 3]],
                               "marker.color": [adata.obs.iloc[:,1]],
                                "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 2",    
                        method="restyle"
                    ),
                    dict(
                        args=[{"x": [adata.obsm['X_diffmap'][:,2]],
                               "y": [adata.obsm['X_diffmap'][:, 3]],
                               "z": [adata.obsm['X_diffmap'][:, 4]],
                               "marker.color": [adata.obs.iloc[:,1]],
                                "marker.colorbar.title": "Pseudotime",
                               "text": [adata.obs.iloc[:,0]]}],
                        label="Components choice 3",    
                        method="restyle"
                    
                    )
                ])
            )
        ]
    )

    # Show the plot
    plot = fig.show(config={"displayModeBar": True}, auto_open=True)
    return(plot)