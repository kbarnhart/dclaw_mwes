"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:17:31 2022

@author: krbarnhart
"""

import os

import cmocean
import dclaw.dplot as cd
import dclaw.dplot as local_dplot
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from dclaw.max_extent import calc_max_extent
from pyclaw.data import Data
from pyclaw.geotools import topotools
from pyclaw.plotters import colormaps, geoplot
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

# import pdb


XLL, YLL, XUR, YUR = calc_max_extent()


import yaml

with open("params.yaml") as file:
    params = yaml.safe_load(file)


def fixup_tick(ax, tick_spacing=1000, fontsize=7):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    xmin = np.ceil(xmin / tick_spacing) * tick_spacing
    ymin = np.ceil(ymin / tick_spacing) * tick_spacing
    xticks = np.arange(xmin, xmax, tick_spacing)
    yticks = np.arange(ymin, ymax, tick_spacing)

    if len(xticks) >= 1:
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        xticklabs = [
            "{:d}$^{{ '{:3d} ' }}$".format(int(n // 1000), int(n % 1000))
            for n in xticks
        ]
        yticklabs = [
            "{:d}$^{{ ' {:3d} ' }}$".format(int(n // 1000), int(n % 1000))
            for n in yticks
        ]

        ax.set_xticklabels(
            xticklabs,
            fontsize=fontsize,
            va="center",
            ha="center",
        )

        ax.set_yticklabels(
            yticklabs,
            fontsize=fontsize,
            rotation="vertical",
            va="center",
            ha="center",
        )


# --------------------------
def setplot(plotdata):
    # --------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    import dclaw.dplot as local_dplot
    from numpy import linspace
    from pyclaw.plotters import colormaps, gaugetools, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # -----------------------------------------
    # Three-Panel Figure
    # -----------------------------------------

    plotfigure = plotdata.new_plotfigure(name="pcolor", figno=0)
    plotfigure.show = True
    plotfigure.kwargs = dict(dpi=300, figsize=(6, 12))

    # ------------------------------
    # Panel 1: Thickness
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="thickness")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]
    plotaxes.title = "Debris Thickness"
    plotaxes.title_with_t = True
    plotaxes.title_t_units = "seconds"
    plotaxes.show = True
    plotaxes.axescmd = "subplot(611)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.depth
    plotitem.imshow_cmap = cmocean.cm.turbid
    plotitem.imshow_norm = colors.LogNorm(0.001, 10, clip=True)
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.gridedges_show = 1
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]
    # ------------------------------
    # Panel 2: Velocity
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="velocity")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]
    plotaxes.title = "Velocity"
    plotaxes.title_with_t = False
    plotaxes.show = True
    plotaxes.axescmd = "subplot(612)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False
    plotitem.colorbar_kwargs = {}

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.velocity_magnitude
    plotitem.imshow_cmap = cmocean.cm.speed
    plotitem.imshow_cmin = 0
    plotitem.imshow_cmax = 10
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.gridedges_show = 1
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]
    
    # ------------------------------
    # Panel 3: static angle
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="b_eroded")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]

    plotaxes.title = "Eroded depth, m"
    plotaxes.title_with_t = False
    plotaxes.show = True
    plotaxes.axescmd = "subplot(613)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.b_eroded
    plotitem.imshow_cmap = cmocean.cm.haline
    plotitem.imshow_cmin = 0
    plotitem.imshow_cmax = 3
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.gridedges_show = 1
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]

  
    # ------------------------------
    # Panel 4 surface slope
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="localslope")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]

    plotaxes.title = "Local slope, degrees"
    plotaxes.title_with_t = False
    plotaxes.show = True
    plotaxes.axescmd = "subplot(614)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.local_slope
    plotitem.imshow_cmap = cmocean.cm.haline_r
    plotitem.imshow_norm = colors.BoundaryNorm([0,1,2,3,4,5,8,11,15], ncolors=256, clip=True)
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.gridedges_show = 1
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]




    # ------------------------------
    # Panel 5: Other
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="sigmae")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]

    plotaxes.title = "(litho-pb)/litho"
    plotaxes.title_with_t = False
    plotaxes.show = True
    plotaxes.axescmd = "subplot(615)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.sigma_e_over_lithostatic
    plotitem.imshow_cmap = cmocean.cm.dense
    plotitem.imshow_norm = colors.LogNorm(0.001, 1, clip=True)    
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.gridedges_show = 1
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]


    # ------------------------------
    # Panel 5: Other
    # ------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name="hchi")
    plotaxes.xlimits = [XLL, XUR]
    plotaxes.ylimits = [YLL, YUR]

    plotaxes.title = "Species 1"
    plotaxes.title_with_t = False
    plotaxes.show = True
    plotaxes.axescmd = "subplot(616)"
    plotaxes.show_gauges = True
    plotaxes.gauge_kwargs = {
        "format_string": "D",
        "markersize": 3,
        "mfc": "white",
        "mec": "black",
        "mew": 0.2,
        "fontsize": 6,
    }

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.plot_var = local_dplot.eta
    plotitem.add_colorbar = False

    # Debris
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = local_dplot.species1_fraction
    plotitem.imshow_cmap = cmocean.cm.oxy
    plotitem.imshow_norm = colors.Normalize(0, 1, clip=True)    
    plotitem.imshow_alpha = 1.0
    plotitem.add_colorbar = True
    plotitem.gridedges_show = 1
    plotitem.colorbar_kwargs = {"fraction": 0.03}
    plotitem.amr_gridedges_show = [False, True, True, True, True]
    plotitem.amr_gridedges_linewidth = [0, 0.2, 0.2, 0.2, 0.2]
    plotitem.amr_gridedges_color = [None, "blue", "red", "green", "magenta"]

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True  # print figures
    plotdata.print_format = "png"  # file format
    plotdata.print_framenos = (
        "all"  # range(70,190,10)  # list of frames to print [0,1,2,3]#
    )
    plotdata.print_gaugenos = "all"  # list of gauges to print
    plotdata.print_fignos = [0]  # list of figures to print
    plotdata.html = False  # create html files of plots?
    plotdata.html_homelink = "../README.html"  # pointer for top of index
    plotdata.latex = False  # create latex file of plots?
    plotdata.latex_figsperline = 2  # layout of plots
    plotdata.latex_framesperline = 1  # layout of plots
    plotdata.latex_makepdf = False  # also run pdflatex?
    plotdata.parallel = True
    plotdata.ffmpeg_movie = True  # make animated mp4 movie with ffmpeg
    plotdata.ffmpeg_name = 'animation'

    return plotdata
