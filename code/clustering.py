import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
import operator
import csv
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform
import plotly.express as px

def get_top100(filepath):
    input_file = os.path.join(filepath, 'result_dir.csv')
    with open(input_file, 'r', newline='') as f_input:
        csv_input = csv.DictReader(f_input)
        data = sorted(csv_input, key=lambda row: (row['Name']))[:100]

    output_file = os.path.join(filepath, 'sorted_result.csv')
    with open(output_file, 'w', newline='') as f_output:    
        csv_output = csv.DictWriter(f_output, fieldnames=csv_input.fieldnames)
        csv_output.writeheader()
        csv_output.writerows(data)
    return output_file


def draw_heatmap(filepath,title_name):

    customer_data = pd.read_csv(filepath)
    data_array = customer_data.iloc[:,1:].values
    labels = customer_data.iloc[:,0].values

    fig = ff.create_dendrogram(data_array, orientation='bottom')

    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Heatmap
    dendro_leaves = fig['layout']['xaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))

    data_dist = pdist(data_array)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]
    heatmap = [
        go.Heatmap(
            x = fig['layout']['xaxis']['tickvals'],
            y = fig['layout']['xaxis']['tickvals'],
            z = heat_data,
            colorscale = 'Blues'
        )
    ]
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width':1000, 'height':1000,
                            'showlegend':False, 'hovermode': 'closest', 'title':title_name
                            })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [0.15, 1],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'ticktext':labels[dendro_leaves],
                                    'tickvals':heatmap[0]['x'],
                                    'ticks':""})
                                    
    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'ticktext':labels[dendro_leaves],
                                    'tickvals':heatmap[0]['x'],
                                    'ticks': ""
                            })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.85, 1],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'showticklabels': False,
                                    'ticks':""})

    # Plot!
    fig.show()

def draw_3dmap(filepath):
    data = pd.read_csv(filepath)
    fig = px.scatter_3d(data, x = 'DIST1', y = 'DIST2', z = 'DIST3', color= 'Name')
    fig.show()


def do_clustering(filepath, num_file):
    original_file = os.path.join(filepath, 'result_dir.csv')
    draw_heatmap(original_file, 'All Crystals')
    draw_3dmap(original_file)
    if num_file > 100:
        top100_file = get_top100(filepath) 
        draw_heatmap(top100_file, 'First 100 crystals')
