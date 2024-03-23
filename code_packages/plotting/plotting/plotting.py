import seaborn as sb
import matplotlib.pyplot as pp
import matplotlib as mp
import numpy as np
import pandas as pd
from typing import List, Tuple

sb.set_context('paper')

CMAPS = {
    'rdbu_paired': [
        (0.6509803921568628, 0.807843137254902, 0.8901960784313725),
        (0.12156862745098039, 0.47058823529411764, 0.7058823529411765),
        (0.984313725490196, 0.6039215686274509, 0.6),
        (0.8901960784313725, 0.10196078431372549, 0.10980392156862745),
    ],
    'green' : sb.light_palette('#125e03', as_cmap=True),
    'purple' : sb.light_palette('#8f0a7d', as_cmap=True),
    'blue_orange': sb.diverging_palette(250, 30, l=65, as_cmap=True),
    'green_purple': sb.diverging_palette(145, 300, s=60, as_cmap=True),
    'ireland': sb.diverging_palette(145, 30, l=65, as_cmap=True),
    'mint_gold': sb.diverging_palette(200, 60, l=65, as_cmap=True),
    'purple_pink': sb.diverging_palette(280, 0, l=40, as_cmap=True),
}

def basic_heatmap(df, annot=True, robust=False, size=None, labels=None, fmt='.2f', bar_range=None, *args, **kwargs):
    sb.set_theme(
        context=kwargs.get('context', 'notebook'),
        style = kwargs.get('style', 'darkgrid'),
        palette = kwargs.get('palette', None),
        font_scale = kwargs.get('font_scale', 1.),
        )
    center = kwargs.get('center', 0.)
    norm = kwargs.get('norm', mp.colors.Normalize())
    cmap = kwargs.get('cmap', 'RdBu')
    xrotation = kwargs.get('xrotation', 0)
    yrotation = kwargs.get('yrotation', 0)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)
    cbar = kwargs.get('cbar', True)
    cbarlabel = kwargs.get('cbarlabel', None)
    annot_kws = kwargs.get('annot_kws', {})
    cbar_kws = kwargs.get('cbar_kws', {})
    if size:
        fig, ax = pp.subplots(figsize=size)
    else:
        fig, ax = pp.subplots()
    if bar_range:
        ax = sb.heatmap(
            df, cmap=cmap, center=center, 
            annot=annot, annot_kws=annot_kws, robust=robust, 
            ax=ax, fmt=fmt, norm=norm, vmin=bar_range[0], vmax=bar_range[1],
            cbar=cbar, cbar_kws=cbar_kws
            )
    else:
        ax = sb.heatmap(
            df, cmap=cmap, center=center, annot=annot, annot_kws=annot_kws, 
            robust=robust, ax=ax, fmt=fmt, norm=norm, cbar=cbar, cbar_kws=cbar_kws)
    if not labels:
        labels = ax.get_xticklabels()
    ax.set_xticklabels(labels=labels, rotation=xrotation)
    ylabels = ax.get_yticklabels()
    ax.set_yticklabels(labels=ylabels, rotation=yrotation)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if cbarlabel:
        ax.collections[0].colorbar.set_label(cbarlabel)
    
    return ax.get_figure(), ax

def basic_clustermap(df, trees=(False, False), ticks=True, figsize=None, cbar_pos=None, labels=True, *args, **kwargs):
    sb.set_theme(
        context=kwargs.get('context', 'notebook'),
        style = kwargs.get('style', 'darkgrid'),
        palette = kwargs.get('palette', None),
        font_scale = kwargs.get('font_scale', 1.),
    )
    col_cluster = kwargs.get('col_cluster', False)
    if figsize:
        ax = sb.clustermap(df, cmap='RdBu', center=0., col_cluster=col_cluster, 
        figsize=figsize, cbar_pos=cbar_pos, xticklabels=labels, fmt='.2f')
    else:
        ax = sb.clustermap(df, cmap='RdBu', center=0., col_cluster=col_cluster, 
        cbar_pos=cbar_pos, xticklabels=labels, fmt='.2f')
    ax.ax_row_dendrogram.set_visible(trees[0])
    ax.ax_col_dendrogram.set_visible(trees[1])
    if not ticks:
        h = ax.ax_heatmap
        h.set_yticklabels([]); h.set_yticks([]);
    return ax

def basic_barplot(x, y, size=None, labels=False, *args, **kwargs):
    palette = kwargs.get('colors', sb.color_palette())
    context = kwargs.get('context', 'notebook')
    style = kwargs.get('style', 'darkgrid')
    font_scale = kwargs.get('font_scale', 1.)
    rotation = kwargs.get('rotation', 0)
    xlabel = kwargs.get('xlabel', None)

    sb.set_theme(
        context = context,
        style = style,
        palette = palette,
        font_scale = font_scale,
    )
    if size:
        fig, ax = pp.subplots(
                figsize=(size[0], size[1])
            )
    else:
        fig, ax = pp.subplots()

    ax.set_xlabel(xlabel)
    if labels:
        ax.set_xticklabels(labels, rotation=rotation);
    else:
        ax.set_xticklabels(x, rotation=rotation)

    sb.barplot(
        x=x, 
        y=y, 
        ax=ax,
        palette=palette,
    )
 
    return fig, ax

def range_map_bar_sections(
    df: pd.DataFrame, 
    df_upper: pd.DataFrame, 
    df_lower: pd.DataFrame, 
    labels=None,
    fontsize=20,
    size=(20,10),
    *args, 
    **kwargs
    ) -> Tuple[mp.figure.Figure, mp.axes.Axes]:
    """
    Plots FVA solutions where the valid flux range is a section of a bar graph for each reaction.
    :param df: DF containing flux values (usually pFBA)
    :type df: pd.DataFrame
    :param df_upper: DF containing upper limit for flux (usually FVA)
    :type df_upper: pd.DataFrame
    :param df_lower: DF containing lower limit for flux (usually FVA)
    :type df_lower: pd.DataFrame
    :param labels: Set of labels to use for legend, defaults to None (inferring from DFs)
    :type labels: List[str] or None, optional
    :param fontsize: Fontsize for axes and legend, defaults to 20
    :type fontsize: int
    :param size: Figure size
    :type size: Tuple[int]
    :raises ValueError: raises ValueError if DF dimensions do not match
    :return: Tuple of matplotlib figure and axis 
    :rtype: Tuple[mp.figure.Figure, mp.axes.Axes]
    kwargs: 
    legend - position of legend (see pyplot.legend) [str]
    colors - color palette (f.ex. from seaborn, or list of hex/rgbi colors)
    rotation - rotation of x-labels in degrees [int]
    width - width of columns [float]
    """
    _check_df_compatibility(df, df_upper, df_lower)

    n_cols = len(df.columns)
    w = kwargs.get('width', 2* (1 / n_cols))
    pos = range(len(df.index))
    colors = kwargs.get('colors', sb.color_palette())
    loc = kwargs.get('legend', 'best')
    rot = kwargs.get('rotation', 0)
    hline = kwargs.get('hline', False)

    # make figure
    fig, ax = pp.subplots(figsize=size)
    if not labels:
        labels = df.columns
    for i, c in enumerate(df.columns):
        x = [y+w*i for y in pos]
        height = df_upper.loc[:, c].values - df_lower.loc[:, c].values + .025
        bottom = df_lower.loc[:, c].values - .0125
        ax.bar(
            x = x, 
            height = height, 
            bottom = bottom,
            label = labels[i],
            width = w,
            color = colors[i],
        )

    # create legend
    ax.legend(
        loc=loc, 
        labels=labels,
        fontsize=fontsize,
    )
    
    # set ticks
    yticks = ax.get_yticks()
    ylabels = [str(tick) for tick in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontdict={'fontsize': fontsize})
    ax.set_xticks(range(len(df.index)))
    ax.set_xticklabels(df.index, rotation=rot, fontdict={'fontsize': fontsize})
    
    # draw lines for y-axis and y-axis intercept
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color('#000000')
    ax.spines['left'].set_linewidth(3)
    if isinstance(hline, float):
        ax.axhline(hline, color='k', linewidth=3)
    
    return fig, ax

def range_map_bar_err(
    df: pd.DataFrame, 
    df_upper: pd.DataFrame, 
    df_lower: pd.DataFrame, 
    size=(20,10), 
    fontsize=20, 
    labels=None, 
    *args, 
    **kwargs
    ) -> Tuple[mp.figure.Figure, mp.axes.Axes]:
    """
    Plots FVA solutions where one solution (usually pFBA) is the bar size and the valid flux range is the error bars.
    :param df: DF containing flux values (usually pFBA)
    :type df: pd.DataFrame
    :param df_upper: DF containing upper limit for flux (usually FVA)
    :type df_upper: pd.DataFrame
    :param df_lower: DF containing lower limit for flux (usually FVA)
    :type df_lower: pd.DataFrame
    :param labels: Set of labels to use for legend, defaults to None (inferring from DFs)
    :type labels: List[str] or None, optional
    :param fontsize: Fontsize for axes and legend, defaults to 20
    :type fontsize: int
    :param size: Figure size
    :type size: Tuple[int]
    :raises ValueError: raises ValueError if DF dimensions do not match
    :return: Tuple of matplotlib figure and axis 
    :rtype: Tuple[mp.figure.Figure, mp.axes.Axes]
    kwargs: 
    legend - position of legend (see pyplot.legend) [str]
    colors - color palette (f.ex. from seaborn, or list of hex/rgbi colors)
    rotation - rotation of x-labels in degrees [int]
    width - width of columns [float]
    hline - where to draw horizontal line indicating zero [float, default: None]
    xlabel - label for the x-axis [str, default: None]
    ylabel - label for the y-axis [str, default: None]
    """
    _check_df_compatibility(df, df_upper, df_lower)
    n_cols = len(df)
    n = range(n_cols)

    width = kwargs.get('width', 2*(1 / n_cols))
    colors = kwargs.get('colors', sb.color_palette())
    loc = kwargs.get('legend', 'best')
    rot = kwargs.get('rotation', 0)
    hline = kwargs.get('hline', None)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)

    fig, ax = pp.subplots(figsize=size)
    if not labels:
        labels = df.columns
    
    # draw plot
    for i, c in enumerate(df.columns):
        vals = df[c].values
        pos = [nx + width*i for nx in n]
        yerr_up = abs(df_upper[c].values - vals)
        yerr_down = abs(df_lower[c].values - vals)
        errs = np.zeros((2, len(df)))
        errs[1, :] = yerr_up; errs[0, :] = yerr_down    
        ax.bar(
            x=pos, 
            height=vals, 
            yerr=errs, 
            label=labels[i], 
            width=width, 
            color=colors[i], 
            ecolor='black',
        )
    # draw legend
    ax.legend(
        loc=loc, 
        labels=labels,
        fontsize=fontsize,
    )
    # set label axes
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fontsize)

    # set ticks
    yticks = ax.get_yticks()
    ylabels = [str(tick) for tick in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontdict={'fontsize': fontsize})
    ax.set_xticks(range(len(df.index)))
    ax.set_xticklabels(df.index, rotation=rot, fontdict={'fontsize': fontsize})

    # draw lines for y-axis and y-axis intercept
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color('#000000')
    ax.spines['left'].set_linewidth(3)
    if isinstance(hline, float):
        ax.axhline(hline, color='k', linewidth=3)

    return fig, ax

def _check_df_compatibility(df, df_upper, df_lower):
    if not len(df.columns) == len(df_upper.columns) == len(df_lower.columns):
        raise ValueError('Column lengths do not match')

def savefig(figure, path):
    for ft in ['.svg', '.eps', '.pdf', '.png']:
        figure.savefig(f'{path}{ft}')

def shorten_index(df: pd.DataFrame, length=20)-> pd.DataFrame:
    tags_short = []
    for tag in df.index:
        if len(tag) > length:
            try:
                tags_short.append(tag[:length] + '[...]')
            except IndexError:
                tags_short.append(tag)
        else: 
            tags_short.append(tag)
    replace = dict(zip(df.index, tags_short))
    
    return df.rename(index=replace)
