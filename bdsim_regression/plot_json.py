import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from math import isqrt
from .plot_common import two_factors
from .plot_common import figure

def compare(json_file1, json_file2, show = False, file_name = "comparison_json.html"):
    json_dict1 = json.load(open(json_file1))
    json_dict2 = json.load(open(json_file2))

    hist_dict1 = json_dict1['histograms']
    hist_dict2 = json_dict2['histograms']

    if len(hist_dict1) != len(hist_dict2):
        raise Exception("Histogram length do not match")

    # number of columns and rows
    nrow, ncol = two_factors(len(hist_dict1))

    # create figure
    fig = figure(nrow, ncol)

    # loop over histograms
    for ihist, hist in enumerate(hist_dict1):
        icol = ihist % ncol
        irow = int((ihist - icol)/nrow)

        icol = icol + 1
        irow = irow + 1
        if hist_dict1[hist]['type'] == "histo1d":
            t1 = hist1d(hist_dict1[hist], hist)
            t2 = hist1d(hist_dict2[hist], hist)
        elif hist_dict1[hist]['type'] == "bar1d" :
            t1 = bar1d(hist_dict1[hist], hist, x_as_string=True)
            t2 = bar1d(hist_dict2[hist], hist, x_as_string=True)

        fig.add_trace(t1, row=irow, col=icol)
        fig.add_trace(t2, row=irow, col=icol)

        fig.update_xaxes(title_text=hist_dict1[hist]['xlabel'], row=irow, col=icol)
        if "yscale" in hist_dict1[hist] :
            fig.update_yaxes(type=hist_dict1[hist]['yscale'], row=irow, col=icol)

    # fig.update_layout(template="plotly_white")
    if show :
        fig.show()

    if file_name is not None  and file_name != "":
        fig.write_html(file_name, include_plotlyjs="cdn")

def bar1d(histDict,
          name = "data",
          x_as_string = False) :

    x = np.array(histDict['xcentres'])
    y = np.array(histDict['contents'])

    if x_as_string :
        x = x.astype(str)

    return go.Bar(
        x=x,
        y=y,
        name=name
    )

def hist1d(histDict,
           name ="data",
           x_as_string = False):

    try :
        x = np.array(histDict['xedges'])
    except :
        x = np.array(histDict['xcentres'])

    if x_as_string :
        x = x.astype(str)

    y = np.append(histDict['contents'], histDict['contents'][-1])

    sca = go.Scatter(
        x=x,
        y=y,  # repeat last value so len(y) == len(edges)
        mode="lines",
        line_shape="hv",  # horizontal-then-vertical steps
        name=name,
    )

    return sca

def hist2d(histDict,
           name = "data") :
    xedges = np.array(histDict['xedges'])
    yedges = np.array(histDict['yedges'])
    contents = np.array(histDict['contents'])

    return go.Heatmap(x=0.5 * (xedges[:-1] + xedges[1:]),  # bin centres
                      y=0.5 * (yedges[:-1] + yedges[1:]),
                      z=contents.T,  # transpose: histogram2d is indexed [x, y], heatmap wants [y, x]
                      # colorscale="Viridis",
                      # name=name
    )
    
    return
