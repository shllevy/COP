from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd
import numpy as np
from sklearn.metrics import auc
import sys
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import scipy.stats as stats
import seaborn as sns
import os
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from sklearn import preprocessing


def sort_antigen_by_dilution(df, antigen, dilution_ind=-2):
    col = df.columns[df.columns.str.startswith(antigen)]
    col_dil = col.map(lambda x: float(x.split('_')[dilution_ind]))
    temp = pd.Series(index = col, data=col_dil).sort_values(ascending=False)
    return temp.index, temp.values



def compute_response_groups_by_antigen_ind (df=pd.DataFrame(), antigen=str, infected_col='sars_cov2_infection', quantiles=[0.25, 0.75], print_flag=False, p_threshold=1):
    """
    Quantile analysis for a single antigen computing response groups - should be moved to amutils (with modifications)
    Booster study
    """
    curr_df = df.copy()
    a = antigen  # Lilach added   
    curr_quantiles = curr_df[a].quantile(quantiles)
    
    low_ptids = curr_df[curr_df[a] <= curr_quantiles[quantiles[0]]].index
    high_ptids = curr_df[curr_df[a] >= curr_quantiles[quantiles[1]]].index

    mid_ptids = curr_df[curr_df[a] > curr_quantiles[quantiles[0]]].index.intersection(curr_df[curr_df[a] < curr_quantiles[quantiles[1]]].index)
    
    col_name = a + '_response_group'
    
    curr_df.loc[low_ptids, col_name] = 'low'
    curr_df.loc[high_ptids, col_name] = 'high'
    curr_df.loc[mid_ptids, col_name] = 'mid'
    
    low_counts = curr_df.loc[low_ptids][infected_col].value_counts()
    high_counts = curr_df.loc[high_ptids][infected_col].value_counts()
    mid_counts = curr_df.loc[mid_ptids][infected_col].value_counts()
    
    #verify all values exist
    for val_count in [low_counts, high_counts, mid_counts]:
        for val in ['yes', 'no']:
            if val in val_count.index:
                continue
            else:
                val_count.loc[val] = 0

    table = [  [high_counts["no"], high_counts["yes"]],[low_counts["no"], low_counts["yes"]] ]
    table_str = "".join(['a=', str(low_counts["no"]), \
                         ', b=', str(low_counts["yes"]), \
                         ', c=', str(high_counts["no"]), \
                         ', d=', str(high_counts["yes"])])
    # print(table_str)
    oddsratio, p_value = scipy.stats.fisher_exact(table)
    
    # dictionary mapping infection rates to response groups:
    infection_rates = {}
    infection_rates['low'] = 100 * low_counts["yes"]/(low_counts["no"] + low_counts["yes"])
    infection_rates['mid'] = 100 * mid_counts["yes"]/(mid_counts["no"] + mid_counts["yes"])
    infection_rates['high'] = 100 * high_counts["yes"]/(high_counts["no"] + high_counts["yes"])
    
    
    if print_flag:
        if (p_value < p_threshold):
            print('----------------------------------------------------------------------------------')
            print('Stats comparing infection rates in {0} low vs. high groups (n = {1}):'.format( a, int(curr_df.shape[0])))
            print('----------------------------------------------------------------------------------')

            print('low')
            print(low_counts)
            print('high')
            print(high_counts)

            print('Fisher test comparing low vs. high infection rates: OR = {0:.3f}, p = {1:.5f}\n'.format(oddsratio, p_value))
            print('Infection rate in low group is: {0:.1f}%'.format(infection_rates['low']))
            print('Infection rate in mid group is: {0:.1f}%'.format(infection_rates['mid']))    
            print('Infection rate in high group is: {0:.1f}%'.format(infection_rates['high']))    
            print('\n')
    return oddsratio, p_value, curr_df, infection_rates

def compute_response_groups_by_antigen(df, isotype, antigen, vaccinated, infected_col='infected', quantiles=[0.25, 0.75], 
                                       print_flag=True, p_threshold=1, time_column='time_point'):
    """
    Quantile analysis for a single antigen computing response groups - should be moved to amutils (with modifications)
    Booster study
    """
    
    a = antigen  # Lilach added   
    curr_df = df[(df.isotype == isotype) & (df[time_column]=='T0') & (df.vaccinated == vaccinated)].copy()
    curr_quantiles = curr_df[a].quantile(quantiles)
    
    low_ptids = curr_df[curr_df[a] <= curr_quantiles[quantiles[0]]].index
    high_ptids = curr_df[curr_df[a] >= curr_quantiles[quantiles[1]]].index

    mid_ptids = curr_df[curr_df[a] > curr_quantiles[quantiles[0]]].index.intersection(curr_df[curr_df[a] < curr_quantiles[quantiles[1]]].index)
    
    col_name = a + '_response_group'
    
    curr_df.loc[low_ptids, col_name] = 'low'
    curr_df.loc[high_ptids, col_name] = 'high'
    curr_df.loc[mid_ptids, col_name] = 'mid'
    
    low_counts = curr_df.loc[low_ptids][infected_col].value_counts()
    high_counts = curr_df.loc[high_ptids][infected_col].value_counts()
    mid_counts = curr_df.loc[mid_ptids][infected_col].value_counts()

    table = [  [high_counts['no'], high_counts['yes']],[low_counts['no'], low_counts['yes']] ]
    table_str = "".join(['a=', str(low_counts['no']), \
                         ', b=', str(low_counts['yes']), \
                         ', c=', str(high_counts['no']), \
                         ', d=', str(high_counts['yes'])])
    #print(table_str)
    oddsratio, p_value = scipy.stats.fisher_exact(table)
    
    # dictionary mapping infection rates to response groups:
    infection_rates = {}
    infection_rates['low'] = 100 * low_counts['yes']/(low_counts['yes'] + low_counts['no'])
    infection_rates['mid'] = 100 * mid_counts['yes']/(mid_counts['yes'] + mid_counts['no'])
    infection_rates['high'] = 100 * high_counts['yes']/(high_counts['yes'] + high_counts['no'])
    
    
    if print_flag:
        if (p_value < p_threshold):
            vac_str = 'vaccinated' if vaccinated == 'yes' else 'un-vaccinated'
            print('Stats comparing {0} infection rates in {1} {2} low vs. high groups (n = {3}):'.format(vac_str, a, isotype, int(curr_df.shape[0])))
            print('----------------------------------------------------------------------------------')

            print('low')
            print(low_counts)
            print('high')
            print(high_counts)

            print('Fisher test comparing low vs. high infection rates: OR = {0:.3f}, p = {1:.5f}\n'.format(oddsratio, p_value))
            print('Infection rate in low group is: {0:.1f}%'.format(infection_rates['low']))
            print('Infection rate in mid group is: {0:.1f}%'.format(infection_rates['mid']))    
            print('Infection rate in high group is: {0:.1f}%'.format(infection_rates['high']))    
            print('----------------------------*-----------------------------------')
            print('\n')
    return oddsratio, p_value, curr_df, infection_rates


def compute_groups_quantiles_for_antigens(df=pd.DataFrame(), isotype=str, antigens=[], quantiles=[0.25, 0.75]):
    """
    Quantile analysis for a single antigen computing response groups - should be moved to amutils (with modifications)
    Booster study

    antigens: list of column names for which quantiles should be calculated
    """

    curr_df = df[df.isotype == isotype].copy()
    for a in antigens:
        for i in curr_df.index:
            curr_df.loc[i, a + '_percentile'] = stats.percentileofscore(curr_df[a], curr_df.loc[i, a],
                                                                        kind='strict') / 100
            # Compute the percentile rank of a score relative to a list of scores.

        curr_quantiles = curr_df[a].quantile(quantiles)

        low_ptids = curr_df[curr_df[a] <= curr_quantiles[quantiles[0]]].index
        high_ptids = curr_df[curr_df[a] >= curr_quantiles[quantiles[1]]].index
        mid_ptids = curr_df[curr_df[a] > curr_quantiles[quantiles[0]]].index.intersection(
            curr_df[curr_df[a] < curr_quantiles[quantiles[1]]].index)

        col_name = a + '_response_group_' + isotype

        curr_df.loc[low_ptids, col_name] = 'low'
        curr_df.loc[high_ptids, col_name] = 'high'
        curr_df.loc[mid_ptids, col_name] = 'mid'

    #     return(curr_df.index, curr_df['percentile'].values/100, curr_df[col_name])
    return curr_df

def compute_breadth_and_magnitude_summary_stats(arr_df, prot_names, prot_strs, ind_dict, pos_threshold=2000, compute_gmean_mag=True, compute_breadth=True):
    """
    Breadth and Magnitude summary statistics for array data:
    compute summary statistics of the array - breadth and magnitude and stores them in
    new columns in arr_df. Uses the ind_dict dictionary for each protein in prot_names.
    also computes the geometric mean of the magnitude (August 2019)

    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    prot_names: List
        list of strings of protein antigens from which peptides are on the array_data_filename
    prot_strs: List
        names of proteins in prot_names used for plotting and vizualization.
    ind_dict: dictionary
        dictionary mapping prot_names to columns in the arr_df
    pos_threshold: int
        threshold to define a positive response for breadth computation.
        Default is set to 2000.
    compute_gmean_mag: Bollean
        do not calculate gmean magnitudes if False.
        default is True.
    compute_breadth: Bollean
        do not calculate breadth if False.
        default is True.
    Returns:
    -------
    arr_df: pandas.DataFrame
        dataframe with magnitude and breadth columns added for each protein in prot_strs.
    """
    for p, s in zip(prot_names, prot_strs):
        # insert new columns into dataframe for overall magnitude for each strain
        arr_df.loc[:, s + '_magnitude'] = arr_df[ind_dict[p]].sum(axis=1)

        if compute_gmean_mag:
            arr_df.loc[:, s + '_gmean_magnitude'] = scipy.stats.gmean(arr_df[ind_dict[p]], axis=1)

        # breadth - response is positive if it is above the mean response of that antigen across all samples
        if compute_breadth:
            for col in ind_dict[p]:
                arr_df.loc[:, col + '_binarized'] = \
                    arr_df[col].map(lambda s: 1 if s > pos_threshold else 0)

            arr_df.loc[:, s + '_breadth'] = \
                arr_df.loc[:, [col + '_binarized' for col in ind_dict[p]]].sum(axis=1).astype(float)

    return arr_df

"""plots"""


def plot_antigen_boxplots_by_infection_stat(df, antigens, antigen_names=None, group_column='group', group_order=None, hue='infected_T1', 
                                            hue_order=None, timepoint='T0', y_scale=None, sharey=True, time_column='time_point', 
                                            figsize=(18,5), showfliers=True, title=None, x_label='', dot_size=2.5, savepath=None, 
                                            fig_prefix=None, colormap=None, annot=None, legend=True, width=0.8):
    """
    Plots specific antigens comparing infected to un-infected at a given timepoint. Each antigen is plotted in a specific subplot.
    The hue column is 'infected'.

    Parameters:
    ----------
    df: pd.DataFrame
        the array dataframe
    antigens: pd.Index | list
        list of antigens to plot (columns in df)
    antigen_names: dictionary | None
        mapping of antigens to antigen_names for plotting. Default is none, in which case the antigens themselves will be used.
    group_column: String
        name of column denoting groups
    group_order: list
        list of groups by specified order. Default is None
    hue_order: list
        list of hue groups by specified order. Default is None
    timepoint: String
        name of timepoint for plotting response of infected vs. un-infected. Default is T0
    y_scale: String:
        if set to None (default) will use original scale. If set to log10, or log2 will use a logarithmic scale for the y-axis
    dot_size: float    
        size of dots in swarmplot - used to make sure all points are plotted when plot is crowded. Default is 2.5
    figsize: tuple
        size of figure
    showfliers: boolean
        If True outliers will be plotted
    savepath: String
        path to save figure. If set to None, will not save figure.
    file_prefix: String
        prefix of filename when saving files

    Returns:
    -------
    f, ax: handles
        figure and axis handles for tweaking figure

    """
    if colormap == None:
        u_cmap = sns.color_palette("Greens", 3)
        i_cmap = sns.color_palette('Reds', 3)
    else:
        u_cmap = sns.color_palette(colormap[0], 3)
        i_cmap = sns.color_palette(colormap[1], 3)
    stats_dict = {0.05:'*', 0.001:'**', 0.0001:'***', 0.00001:'****'}
    xtick_dict = {'no':'  Three \n doses', 'yes':' Four \n doses'}
    if hue_order == None:
        h_order = df[hue].unique()
    elif hue_order != None:
        h_order = hue_order

    if group_order == None:
        g_order = df[group_column].unique()
    elif group_order != None:
        g_order = group_order
    f, ax = plt.subplots(1,len(antigens), figsize=figsize,sharey=sharey)

    # multiply values by 100 to get positive values on log scale:
    if (y_scale == 'log10') or (y_scale == 'log2'):
        df.loc[:, antigens] = df[antigens] * 1000
        
    if antigen_names is None:
        antigen_names = antigens

    for ax_i, a in enumerate(antigens):
        curr_ax = ax.flatten()[ax_i]
        if y_scale == 'log10':
            sns.swarmplot(data=df[df[time_column] == timepoint], x=group_column, y=np.log10(df[a]), hue=hue, order=g_order, hue_order=h_order,  ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=False,
                medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 2},
                whiskerprops={'visible': False},
                zorder=10,
                data=df, x=group_column, y =np.log10(df[a]), hue=hue, order=g_order, hue_order=h_order,
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=curr_ax, width=width)
        elif y_scale == 'log2':
            sns.swarmplot(data=df[df[time_column] == timepoint], x=group_column, y=np.log2(df[a]), hue=hue, order=g_order, hue_order=h_order, ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=False,
                medianprops={'visible': True, 'color': 'k', 'ls': '-'},
                whiskerprops={'visible': True},
                zorder=10,
                data=df, x=group_column, y=np.log2(df[a]), hue=hue,  order=g_order, hue_order=h_order,
                showfliers=False,
                showbox=False,
                showcaps=True,
                ax=curr_ax)
        else:
            sns.swarmplot(data=df[df[time_column] == timepoint], x=group_column, y=df[a], hue=hue, order=g_order, hue_order=h_order,  ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=False,
                medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 2},
                whiskerprops={'visible': False},
                zorder=10,
                data=df, x=group_column, y=df[a], hue=hue,  order=g_order, hue_order=h_order, 
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=curr_ax)
            
        if annot == True:              
            for g, g_pos in zip(g_order,np.arange(0,1,1/len(g_order))):
                st, pval = scipy.stats.ranksums(df[(df[time_column]==timepoint)&(df[group_column]==g)&(df[hue] == h_order[0])][a], 
                                                df[(df[time_column]==timepoint)&(df[group_column]==g)&(df[hue] == h_order[1])][a])
                for p, s_pos in zip([0.00001, 0.0001, 0.001, 0.05], [0.162,0.184,0.222,0.228]):
                    if pval < p:
                        curr_ax.text(g_pos+s_pos, 0.9, stats_dict[p], fontsize=17,transform=curr_ax.transAxes)
                        curr_ax.hlines(0.9, g_pos+0.15, g_pos+0.35, colors='black', linestyles='-',transform=curr_ax.transAxes)
                        break

#         curr_ax.set_ylim(3.6, 14.5)
       
        curr_ax.set_ylabel('')
        curr_ax.set_xlabel('')
        curr_ax.set_xticklabels([xtick_dict[g_order[0]], xtick_dict[g_order[1]]])
        curr_ax.set_title(antigen_names[a])
        
        # how to get the collection of Artists by which we can change the color of each box in our boxplot
        c = curr_ax.collections
        
        c[0].set_facecolor(u_cmap[2])
        c[1].set_facecolor(i_cmap[2])
        c[2].set_facecolor(u_cmap[2])
        c[3].set_facecolor(i_cmap[2])
        curr_ax.get_legend().remove()

    # collect the labels of all lines and boxes in our plot
    lines_labels = [ax.get_legend_handles_labels() for ax in f.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    # change colors of some of the specific handles
    lines[0].set_facecolor(u_cmap[2])
    lines[1].set_facecolor(i_cmap[2])

    lines[0].set_edgecolor(u_cmap[2])
    lines[1].set_edgecolor(i_cmap[2])
    ax[0].set_ylabel('Magnitude MFI')
    # if isotype == 'IgG':
    #     ax[0].set_ylim([8, 12])
    #     ax[1].set_ylim([8, 12])
        
    
    # Using legend on the fig instead of axis generates one legend for all graphs!
    # finally we invoke the legend with the specific lines we want to present - notice location is outside 
    
    if legend == True:
        f.legend([lines[0], lines[1]], ['Uninfected', 'Infected'],  loc = 2, bbox_to_anchor = (1,0.83))
    if title != None: 
        f.suptitle(title, x= 0.527, y=0.88)
    plt.tight_layout()
    if savepath != None:
        plt.savefig(os.path.join(savepath, fig_prefix + '_by_vax_group_boxplots.png'), dpi=300, bbox_inches='tight', facecolor='w')
    return f, ax


def plot_antigen_boxplots_by_timepoints(df, antigens, antigen_names=None, group_column='group', group_order=None,
                                        hue_order=None,
                                        time_column='time_point', y_scale=None, sharey='True', figsize=(18, 5),
                                        annot=None, legend=True, dot_size=4,
                                        showfliers=True, title=None, x_label='', savepath=None, fig_prefix=None,
                                        colormap=None, fig_format='png'):
    """
    Plots specific antigens comparing groups across time. Each antigen is plotted in a specific subplot.

    Parameters:
    ----------
    df: pd.DataFrame
        the array dataframe
    antigens: pd.Index | list
        list of antigens to plot (columns in df)
    antigen_names: dictionary | None
        mapping of antigens to antigen_names for plotting. Default is none, in which case the antigens themselves will be used.
    group_column: String
        name of column denoting groups
    group_order: list
        list of groups by specified order. Default is None
    time_colum: String
        name of time column for plotting response of each group across time
    y_scale: String:
        if set to None (default) will use original scale. If set to log10, or log2 will use a logarithmic scale for the y-axis
    figsize: tuple
        size of figure
    showfliers: boolean
        If True outliers will be plotted
    savepath: String
        path to save figure. If set to None, will not save figure.
    file_prefix: String
        prefix of filename when saving files

    Returns:
    -------
    f, ax: handles
        figure and axis handles for tweaking figure

    """
    stats_dict = {0.05: '*', 0.001: '**', 0.0001: '***', 0.00001: '****'}

    if colormap == None:
        d0_cmap = sns.color_palette("deep")
        d30_cmap = sns.color_palette('bright')
    else:
        d0_cmap = colormap['d0_cmap']
        d30_cmap = colormap['d30_cmap']
    xtick_dict = {'no': '  Three \n doses', 'yes': ' Four \n doses'}

    f, ax = plt.subplots(1, len(antigens), figsize=figsize, sharey=sharey)
    plt.subplots_adjust(wspace=0, hspace=None)
    if antigen_names is None:
        antigen_names = antigens

    for ax_i, a in enumerate(antigens):
        curr_ax = ax.flatten()[ax_i]
        if y_scale == 'log10':
            sns.swarmplot(data=df, x=group_column, y=np.log10(df[a]), hue=time_column, order=group_order,
                          hue_order=hue_order, ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=True,
                        medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 2},
                        whiskerprops={'visible': False},
                        zorder=10,
                        data=df, x=group_column, y=np.log10(df[a]), hue=time_column, order=group_order,
                        hue_order=hue_order,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=curr_ax)
        elif y_scale == 'log2':
            sns.swarmplot(data=df, x=group_column, y=np.log2(df[a]), hue=time_column, order=group_order,
                          hue_order=hue_order, ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=True,
                        medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 2},
                        whiskerprops={'visible': False},
                        zorder=10,
                        data=df, x=group_column, y=np.log2(df[a]), hue=time_column, order=group_order,
                        hue_order=hue_order,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=curr_ax)
        else:
            sns.swarmplot(data=df, x=group_column, y=df[a], hue=time_column, order=group_order, hue_order=hue_order,
                          ax=curr_ax, dodge=True, size=dot_size)
            # plot the median line
            sns.boxplot(showmeans=False,
                        medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 2},
                        whiskerprops={'visible': False},
                        zorder=10,
                        data=df, x=group_column, y=df[a], hue=time_column, order=group_order, hue_order=hue_order,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=curr_ax)

        if hue_order == None:
            h_order = df[time_column].unique()
        elif hue_order != None:
            h_order = hue_order

        if group_order == None:
            g_order = df[group_column].unique()
        elif group_order != None:
            g_order = group_order

        if annot == True:
            
            for g, g_pos in zip(g_order,np.arange(0,1,1/len(g_order))):
                st, pval = scipy.stats.ranksums(df[(df[time_column]==h_order[0])&(df[group_column]==g)][a], 
                                                df[(df[time_column]==h_order[1])&(df[group_column]==g)][a])
                for p, s_pos in zip([0.00001, 0.0001, 0.001, 0.05], [0.18,0.201,0.222,0.228]):
                    if pval < p:
                        curr_ax.text(g_pos+s_pos, 0.9, stats_dict[p], fontsize=17,transform=curr_ax.transAxes)
                        curr_ax.hlines(0.9, g_pos+0.15, g_pos+0.35, colors='black', linestyles='-',transform=curr_ax.transAxes)
                        break


        curr_ax.set_ylabel('')
        curr_ax.set_xlabel(x_label)
        curr_ax.set_xticklabels([xtick_dict[g_order[0]], xtick_dict[g_order[1]]])
        #         curr_ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        curr_ax.set_title(antigen_names[a])

        # how to get the collection of Artists by which we can change the color of each box in our boxplot
        c = curr_ax.collections

        c[2].set_facecolor(d0_cmap[0])
        c[3].set_facecolor(d30_cmap[0])
        c[0].set_facecolor(d0_cmap[1])
        c[1].set_facecolor(d30_cmap[1])
        curr_ax.get_legend().remove()

    ax[0].set_ylabel('Magnitude [MFI]')
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    # collect the labels of all lines and boxes in our plot
    lines_labels = [ax.get_legend_handles_labels() for ax in f.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    # change colors of some of the specific handles
    lines[0].set_facecolor(d0_cmap[0])
    lines[1].set_facecolor(d30_cmap[0])
    lines[4].set_facecolor(d0_cmap[1])
    lines[5].set_facecolor(d30_cmap[1])

    lines[0].set_edgecolor(d0_cmap[0])
    lines[1].set_edgecolor(d30_cmap[0])
    lines[4].set_edgecolor(d0_cmap[1])
    lines[5].set_edgecolor(d30_cmap[1])

    # Using legend on the fig instead of axis generates one legend for all graphs!
    # finally we invoke the legend with the specific lines we want to present - notice location is outside

    if legend == True:
        f.legend([lines[4], lines[5], lines[0], lines[1]],
                 ['Three doses day 0', 'Three doses day 30', 'Four doses day 0', 'Four doses day 30'], loc=2,
                 bbox_to_anchor=(1, 0.83))
    if title != None:
        f.suptitle(title, x=0.55, y=0.88)
    plt.tight_layout()
    if savepath != None:
        plt.savefig(os.path.join(savepath, fig_prefix + '_group_time_boxplots.{0}'.format(fig_format)), dpi=300,
                    bbox_inches='tight', facecolor='w')
    return f, ax

def show_values_on_bars(axs):
    def _show_on_single_plot(ax):
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height() + 2
            value = '{:.1f}'.format(p.get_height())
            ax.text(_x, _y, str(value) + '%', ha="center", size=16)

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)


def plot_sorted_percentile(df, column, colors=['#438380', '#A9C967', '#A8383B'], div=25, figsize=(13, 8)):
    temp_df = df.copy()
    for i in temp_df.index:
        temp_df.loc[i, column + '_percentile'] = stats.percentileofscore(temp_df[column], temp_df.loc[i, column],
                                                                         kind='strict')

    temp_df.loc[temp_df[temp_df[column + '_percentile'] <= div].index, 'color'] = colors[-1]
    temp_df.loc[temp_df[temp_df[column + '_percentile'] >= 100 - div].index, 'color'] = colors[0]
    temp_df['color'] = temp_df.color.fillna(colors[1])
    temp_df = temp_df.sort_values(column + '_percentile', ascending=False)
    color_list = temp_df.color.tolist()

    f, ax = plt.subplots(1, 1, figsize=figsize)
    temp_df[column].plot.bar(color=color_list)
    ax.set_xticklabels([])
    ax.set_title('SARS_CoV_2 Variants magnitude')
    ax.set_ylabel('Magnitude')
    return f, ax


def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarAxes(PolarAxes):

        name = 'radar'

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("unknown value for 'frame': %s" % frame)

        def draw(self, renderer):
            """ Draw. If frame is polygon, make gridlines polygon-shaped """
            if frame == 'polygon':
                gridlines = self.yaxis.get_gridlines()
                for gl in gridlines:
                    gl.get_path()._interpolation_steps = num_vars
            super().draw(renderer)


        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)


                return {'polar': spine}
            else:
                raise ValueError("unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


def collect_spider_data_LF(arr_df, curr_ptids, curr_columns):
    # function for collecting spider data for a given set of ptids, modified by Lilach
    
    # standardize data - to normalize across strains:
    curr_df = arr_df[curr_columns]
    min_max_scaler = preprocessing.MinMaxScaler()
    curr_data_scaled = min_max_scaler.fit_transform(curr_df[curr_columns].values)
    spider_data = pd.DataFrame(index=curr_df.index, columns=curr_df.columns, data=curr_data_scaled)
    title_strs = []
    for p in curr_ptids:

        # generate plot title (including ptid and infection status)
        
        curr_ptid_inds = arr_df.index[arr_df.index.str.startswith(p)]
        title = p
        title_strs.append(title)
    return spider_data, title_strs


def plot_spiderplot_specific_colors(data, curr_timepoints, all_timepoints, 
                    spoke_labels=None, title=None, legend=None, rgrids=[0.2, 0.4, 0.6, 0.8], 
                    donut_sizes=None, donut_labels=None, donut_colors=None, edgecolor=None,
                    spoke_annotation_inds= None, figsize=(6, 6), figpath=None, figprefix=None, start_angle=90,
                    ptid_colors_dict=None, alpha=0.25, line_width=1, max_value=1, chart_radius=1.95, title_height = 1.2, diff_circles=0.1,
                    title_font_size=20, ax_handle=None, donut_labels_textprop=None):
    """
    Plots spider plots for the given list of ptids (can be of length 1) using the spiderplot.py library
    Also for cases when the number of timepoints is not identical for all samples
    (edited by LF)
    
    Parameters:
    ----------
    data: [list | array | series]
        values of the responses to each of the antigens on the web. All timepoints in the same data
    timepoints: [list]
        curr_timepoints names of the data of current ptid (e.g. ['pre', 'post'])
    all_timepoints: [list]
        all timepoint names in the experiment (e.g. ['pre', 'post'])
    spoke_labels: [list | None]
        list of strings for the labels of each of the strains, default is None.
    title: [string | None]
        title fo the plot. Default set to None
    rgrids: [list]
        list of grid values on the web. Default set to [0.2, 0.4, 0.6, 0.8]
    spoke_annotation_inds: [list]
        positions of specidic spokes to label in spoke_labels
    ptid_colors_dict:  [dict]
        e.g. {'pre': sbn.color_palette("tab10")[0], 'post': sbn.color_palette("tab10")[1]}
    alpha: [float 0-1]
        transpancy of the facecolor (inner color of the spider)
    max_value:   [int]
        For normalized spiders =  1.  For un-normalized results = max of spider data values. 
        Reguarly, normalized values are in the range 0-1.
        If the the max value is different, spoke_annotation_position (the positions to get the dot on the spokes circle around the spiders) and rgrids should be updated, 
        so the max value should be specified. 
    chart_radius:   [float]
        Radius of the outer circle
    title_height:   [float]
        The height of the title aboce the chart. A bigger number will give a higher title. In most cases, should be between 1-1.3
    title_font_size; [int]
        Font size of the title
    diff_circles:   [float]
        The space between the inner (radar) and outer (donut) circles
    donut_labels_textprop :[dict | None]
        default is None , dict for text properties - Dict of arguments to pass to the text objects
        https://matplotlib.org/stable/tutorials/text/text_props.html
        :arg 'size' control donut labels size

    Returns:
    -------
    fig: figure
        handle to figure.

    """
    N = data.shape[1]
    theta = radar_factory(N, frame='polygon')

    spoke_annotation_position = max_value+0.05*max_value
    if max_value != 1:
        rgrids = [max_value*x for x in rgrids]


    if ptid_colors_dict is None:
        color_number = all_timepoints.shape[0]
        ptid_colors_dict = dict(zip(all_timepoints, sbn.color_palette("tab10")[0:color_number]))
        

    if ax_handle is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='radar', zorder=2))
        fig.subplots_adjust(top=0.85, bottom=0.05)
    else:
        ax = ax_handle
        fig = ax.get_figure()
        #fig.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.7, wspace=0.5)


    ax.set_title(title,  position=(0.5, 1.25), ha='center', y = title_height, fontsize=title_font_size)

    # plot actual data
    for d, t in zip(data, curr_timepoints):
        c = ptid_colors_dict[t]
        if edgecolor:
            line = ax.plot(theta, d, color=edgecolor, lw=line_width)
        else:
            line = ax.plot(theta, d, color=c, lw=line_width)
        ax.fill(theta, d,  alpha=alpha, color=c)

    if spoke_labels:
        ax.set_varlabels(spoke_labels)

    ax.set_rgrids(rgrids)
    ax.set_yticklabels([])


    if legend is not None:
        ax.legend(legend)
 
    # If provided will plot donut plot around our spider plot for labels of subtypes:
    if donut_sizes is not None:
        if ax_handle is not None:
            ax2 = ax.inset_ax = inset_axes(ax, width='100%', height='100%', loc=10)
            ax2.set_zorder(0)
        else:
            ax2 = fig.add_subplot(111, label="pie axes", zorder=0)
        d = ax2.pie(donut_sizes, labels=donut_labels, colors=donut_colors, radius=chart_radius, autopct=None, startangle=start_angle, labeldistance=1.05,
                    textprops=donut_labels_textprop)
        ax2.set(aspect="equal")

        # add padded white circle to separte donut from spider plot
        if ax_handle is not None:
            ax3 = ax.inset_ax = inset_axes(ax, width='100%', height='100%', loc=10)
            ax3.set_zorder(1)
        else:
            ax3 = fig.add_subplot(111, label='pie', zorder=1)
        d = ax3.pie([1], labels=[''], radius=chart_radius-diff_circles, autopct=None, startangle=start_angle, labeldistance=1.05, colors=['white'])
        ax3.set(aspect="equal")

        ax.set_zorder(2)
    
    # add annotation so specific spokes on the web (used for marking things such as vaccine strains etc.)
    # assumes donut_colors are not None!
    if spoke_annotation_inds is not None:
        
        r = np.arange(0,1, 1/N)
        theta = 2*np.pi * r
        for i, ind in enumerate(spoke_annotation_inds):
            thisr, thistheta = spoke_annotation_position, theta[ind]
            ax.plot([thistheta], [thisr], 'o', color=donut_colors[i])

    if figpath is not None:
        filename = figpath + '/' + figprefix + '_spiderplot'
        fig.savefig(filename)

    return fig


def data_transform(df=pd.DataFrame(), data_value_x=str , data_value_y=list, hue_column=None, ptid_recognition=str):
    """
    transform data to set all columns variable to one column by considering other necessary data
    variable

    :param df: DataFrame containing all data
    :param data_value_x: "str" column including x axis data
    :param data_value_y: "list" of columns including y axis variables
    :param hue_column: "str" column including hue data
    :param ptid_recognition: "str" column including the repetitive information that right for all other conditions
    :return: pandas DataFrame
    """
    if hue_column != None:
        columns = pd.Index(['old_column',ptid_recognition, hue_column, 'variable', data_value_x])

    else:
        columns = pd.Index(['old_column', ptid_recognition,  'variable', data_value_x])

    new_df = pd.DataFrame(columns=columns)

    for col in data_value_y:
        temp_df = df.loc[:, columns.drop(['old_column','variable']).append(pd.Index([col]))]
        temp_df.loc[:, 'old_column'] = col
        temp_df.loc[:, 'x_data'] = temp_df.loc[:, 'old_column'] + '_' + temp_df.loc[:, data_value_x]
        new_df = pd.concat([new_df, temp_df.rename(columns={col:'variable'})])

    new_df = new_df.reset_index()
    return new_df


def get_x_coordinates_of_seaborn_boxplot(ax):
    """

    :param ax: matplotlib.axes.Axes
    :return: x tick positions
    """
    x_coordinates = []
    for c in ax.get_children():
        if type(c) == matplotlib.lines.Line2D:
            if len (c.get_xdata()) == 2:
                x_coordinates = x_coordinates+[((c.get_xdata()[0] + c.get_xdata()[1])/2)]
    return pd.Index(x_coordinates).unique()


def find_repetitive_without_hue(data=None, x_col=None, y_cols=None, recognition=None):
    """

    :param data: DataFrame containing data
    :param x_col: x column including x axis data
    :param y_cols: list, y cols including y axis data
    :param recognition: "str" column including the repetitive information that right for all other conditions
    :return: dictionary including repetitive recognition (set as values) for different conditions (keys)
    """
    final_dict = {}
    set_dict = {}
    set_list = []

    if len(y_cols) == 1:
        for i_, x_ in enumerate(data[x_col].unique()):
            set_dict['{0},{1}'.format(y_cols[0], str(x_))] = set(data[data[x_col] == x_][recognition].values)
            set_list.append('{0},{1}'.format(str(y_cols[0]), str(x_)))

        for set_pos in range(len(set_list) - 1):
            curr_set = set_list[set_pos]
            final_dict[curr_set] = set_dict[curr_set] & set_dict[set_list[set_pos + 1]]
        final_dict[set_list[-1]] = set_dict[set_list[-1]] & set_dict[set_list[-2]]
    else:

        for i_, x_ in enumerate(data[x_col].unique()):
            temp_list = []
            for col in y_cols:
                set_dict['{0},{1}'.format(str(col), str(x_))] = set(
                    data[(data[x_col] == x_) & (data['old_column'] == col)][recognition].values)
                temp_list.append('{0},{1}'.format(str(col), str(x_)))
            set_list.append(temp_list)

        for set_pair in set_list:
            for set_pos in range(len(set_pair) - 1):
                curr_set = set_pair[set_pos]
                final_dict[curr_set] = set_dict[curr_set] & set_dict[set_pair[set_pos + 1]]
            final_dict[set_pair[-1]] = set_dict[set_pair[-1]] & set_dict[set_pair[-2]]

    return (final_dict)


def find_repetitive_with_hue(data=None, y_cols=None, x_col=None, hue_col=None, h_var=None, recognition=None):
    """

    :param data: DataFrame containing data
    :param y_cols: list, y cols including y axis data
    :param x_col: x column including x axis data
    :param hue_col: column including hue data
    :param h_var: order of hue
    :param recognition: "str" column including the repetitive information that right for all other conditions
    :return: dictionary including repetitive recognition (set as values) for different conditions (keys)
    """
    final_dict = {}
    set_dict = {}
    set_list = []

    for i_, x_ in enumerate(data[x_col].unique()):
        for col in y_cols:
            temp_list = []
            for h_ in h_var:
                set_dict['{0},{1},{2}'.format(str(col), str(x_), str(h_))] = set(data[(data[x_col] == x_) &
                                                                                      (data['old_column'] == col) & (
                                                                                              data[hue_col] == h_)][
                                                                                     recognition].values)
                temp_list.append('{0},{1},{2}'.format(str(col), str(x_), str(h_)))
            set_list.append(temp_list)

    for set_ in set_list:
        for set_pos in range(len(set_) - 1):
            final_dict[set_[set_pos]] = set_dict[set_[set_pos]] & set_dict[set_[set_pos + 1]]
        final_dict[set_[-1]] = set_dict[set_[-1]] & set_dict[set_[-2]]
    return (final_dict)


def create_pair_list(x_order=None, y=None, hue=None, hue_order=None):
    """

    :param x_order: list, order of x axis data
    :param y: list, order of y axis data as iy shown on x axis
    :param hue: str, column including hue data
    :param hue_order: order of hue data as it plot on x axis
    :return: array including lists of paired axis
    """

    pair_list = []
    for order in x_order:
        if hue==None:
            for col_i in np.arange(len(y)-1):
                pair_list.append(['{0},{1}'.format(y[col_i], order), '{0},{1}'.format(y[col_i+1], order)])
        else:
            for col in y:
                for hue_i in np.arange(len(hue_order)-1):
                    pair_list.append(['{0},{1},{2}'.format(col, order, hue_order[hue_i]),
                                      '{0},{1},{2}'.format(col, order, hue_order[hue_i+1])])
    return(pair_list)


def create_color_dict(data=None, x_col=None, colors=None, box_list=None, index_dict=None):
    """

    :param data: DataFrame including all data
    :param x_col: str, column including x axis data
    :param colors: str, list of colors as number of x paramters, or dictionary with colors as keys and index as value
    :param box_list: list, order of data in x axis
    :return: dictionary of index as key and color as value
    """
    if type(colors) == list:
        dot_color_dict = {c: index_dict[c_l] for c, c_l in zip(colors, box_list)}
    elif type(colors) == str:
        dot_color_dict = {colors: data.index}
    elif type(colors) == dict:
        dot_color_dict = {k_c: data[data[x_col].isin(colors[k_c])].index for k_c in colors.keys()}

    r_color_dict = {}
    for color, i in dot_color_dict.items():
        r_color_dict.update({i_: color for i_ in i})

    return (r_color_dict)


def pair_boxplot(data=None,
                 y=None,
                 x=None,
                 order=None,
                 y_scale=None,
                 hue=None,
                 hue_order=None,
                 backbone='boxplot',
                 width=0.3,
                 ptid_recognition='ptid',
                 figsize=(33, 14),
                 jitter_scale=None,
                 c_plaette='black',
                 alpha=0.4,
                 ax=None,
                 line_color='black',
                 pair_line_style='--',
                 pair_line_width=0.5,
                 pair_line_alpha=1,
                 median_color='red',
                 median_style='-',
                 median_width=2,
                 pair_flag=True):
    """

    :param data: DataFrame including all data
    :param y: list, of y columns
    :param x: str, x axis data
    :param order: list, order of x axis data
    :param y_scale: None(default) , 'log2'/ 'log10'
    :param hue: str, hue column name (None(default))
    :param hue_order: list of hue order
    :param backbone: str, boxplot / median / None
    :param width: boxplot width 0.9 as default
    :param ptid_recognition: "str" column including the repetitive information that right for all other conditions
    :param figsize: (x axis size, y axis size)
    :param jitter_scale: float, jitter scale of scatter plot
    :param c_plaette: str, list of colors as number of x paramters, or dictionary with colors as keys and index as value
    :param alpha: float
    :param ax: None ax default, ax in case of subplots
    :param line_color: str or dictionary
    :param pair_line_style: '--' as default
    :param pair_line_width: float
    :param pair_line_alpha: float
    :param median_color: color of median line
    :param median_style: '_' as default
    :param median_width: float, 2 as default
    :param pair_flag: True as default. False for plot without lines
    :return: matplotlib.axes.Axes
    """

    if ax == None:
        f, ax = plt.subplots(1, figsize=figsize)
 
    # set defualts:
    if type(line_color) == dict:
        line_color_dict = dict()
        for k in line_color.keys():
            for v in line_color[k]:
                line_color_dict[v] = k
    else:
        line_color_dict = dict(zip(data.index, [line_color] * len(data.index)))

    # test function parameters, for multipule columns as y perform data transformation

    new_df = data_transform(data, x, y, hue, ptid_recognition)
    new_df['variable'] = pd.to_numeric(new_df['variable'])
    new_df.rename(columns={'index': 'old_index'}, inplace=True)
    # display(new_df)

    if y_scale != None:
        new_df.loc[:, 'variable'] = y_scale(new_df.loc[:, 'variable'])

    # set color plate
    if c_plaette == None:
        c_plaette = 'gray'
    else:
        c_plaette = c_plaette

    # set jitter
    if jitter_scale == None:
        jitter_scale = 0.03
    else:
        jitter_scale = jitter_scale

    if order == None:
        order_box = data[x].unique()
    else:
        order_box = order

    order_list_box = []
    final_order_list_box = []

    for o_box in order_box:
        for y_box in y:
            order_list_box.append(','.join([y_box, str(o_box)]))

    if hue == None:
        hue_order = None
        final_order_list_box = order_list_box
        repetitive_dict = find_repetitive_without_hue(data=new_df, y_cols=y, x_col=x, recognition=ptid_recognition)

    elif hue != None:
        if hue_order == None:
            hue_order = data[hue].unique()
        else:
            hue_order = hue_order

        repetitive_dict = find_repetitive_with_hue(data=new_df, y_cols=y, x_col=x,
                                                   recognition=ptid_recognition, hue_col=hue, h_var=hue_order)
        # display(repetitive_dict)
        for x_o in order_list_box:
            for h_order in hue_order:
                final_order_list_box.append('{0},{1}'.format(x_o, h_order))

    if len(y) > 1:
        pair_list = create_pair_list(x_order=order_box, y=y, hue=hue, hue_order=hue_order)
    else:
        pair_list = create_pair_list(x_order=y, y=order_box, hue=hue, hue_order=hue_order)
        new_pair_list = []
        for p in pair_list:
            temp = []
            for p_v in p:
                temp.append((',').join(p_v.split(',')[::-1]))
            new_pair_list.append(temp)
        pair_list = new_pair_list

    # draw boxplot

    if backbone == 'boxplot':
        box = True
        caps = True
        medianprops = {'visible': True, 'color': median_color, 'ls': median_style, 'lw': median_width}
        whiskerprops = {'visible': True}

    elif backbone == 'median':
        box = False
        caps = False
        medianprops = {'visible': True, 'color': median_color, 'ls': median_style, 'lw': median_width}
        whiskerprops = {'visible': False}

    elif backbone == None:
        box = False
        caps = False
        medianprops = {'visible': False}
        whiskerprops = {'visible': False}

  
    sns.boxplot(data=new_df, y=new_df['variable'].values, order=[s.replace(',', '_') for s in order_list_box], hue=hue, x='x_data',
                hue_order=hue_order, color='w', width=width,
                zorder=1, showfliers=False, showbox=box, showcaps=caps, medianprops=medianprops,
                whiskerprops=whiskerprops, ax=ax)
    for patch in ax.artists:
        patch.set_facecolor('None')
    #     patch.set_edgecolor('black')

    coords = get_x_coordinates_of_seaborn_boxplot(ax)

    pos_dict = dict(zip(final_order_list_box, coords))
    index_dict = {}
    jitter_dict = {}
    color_list = []

    for box in final_order_list_box:
        curr_box = box.split(',')
        if len(curr_box) == 2:
            # index_dict[box] = new_df[(new_df.old_column == curr_box[0]) & (new_df[x] == curr_box[1]) & (
            #     new_df[ptid_recognition].isin(repetitive_dict[box]))].index
            index_dict[box] = new_df[(new_df.old_column == curr_box[0])& (new_df[x] == curr_box[1])].index
        elif len(curr_box) == 3:
            # index_dict[box] = new_df[(new_df.old_column == curr_box[0]) & (new_df[x] == curr_box[1]) &
            #                          (new_df[hue] == curr_box[2]) & (
            #                              new_df[ptid_recognition].isin(repetitive_dict[box]))].index
            index_dict[box] = new_df[(new_df.old_column == curr_box[0]) & (new_df[x] == curr_box[1]) &
                                     (new_df[hue] == curr_box[2])].index
        jitter_dict[box] = np.random.normal(loc=pos_dict[box], scale=jitter_scale,
                                            size=len(index_dict[box]))  # set jitter
        if str(box) not in color_list:
            color_list.append(box)

    r_color_dict = create_color_dict(data=new_df, x_col='old_index', colors=c_plaette, box_list=color_list, index_dict=index_dict)
    new_df['dot_color'] = new_df.index.map(r_color_dict).fillna('gray')

    for key in index_dict.keys():
        dot_colors = new_df.loc[index_dict[key], 'dot_color'].values.tolist()
        label = ('\n ').join(key.split(','))
        for i_key in np.arange(len(index_dict[key])):
            ax.plot(jitter_dict[key][i_key], new_df.loc[index_dict[key][i_key], 'variable'], 'o', alpha=alpha,
                    zorder=2, ms=8, mew=1,
                    color=dot_colors[i_key], label=label)
            
    pos_df = pd.DataFrame(index = new_df.index, columns=('x', 'jitter', 'y', 'ptid', 'old_index'))
    for k in index_dict.keys():
        for v_i, v_j in zip(index_dict[k], jitter_dict[k]):
            pos_df.loc[v_i, ['x', 'jitter']] =[k, v_j]
            pos_df.loc[v_i, ['y', 'ptid', 'old_index']] = new_df.loc[v_i, ['variable', 'ptid', 'old_index']].values
            pos_df['x'] = pd.Categorical(pos_df['x'],categories=final_order_list_box,ordered=True)


    if pair_flag != False:
        for p_id in pos_df.ptid.unique():
            curr_df = pos_df[pos_df.ptid == p_id].copy().sort_values('x')
            if len(curr_df) > 2:
                pos_l = np.arange(len(curr_df)-1)
            # elif len(curr_df) == 2:
            #     pos_l = [0]
            else: 
                print(p_id, ' shown only once')
                pos_l = []
            for i_l in pos_l:
                color_idx = new_df.iloc[i_l]['old_index']
                line_c = line_color_dict[color_idx] if color_idx in line_color_dict.keys() else 'black'
                ax.plot([curr_df.iloc[i_l]['jitter'], curr_df.iloc[i_l+1]['jitter']], [curr_df.iloc[i_l]['y'], curr_df.iloc[i_l+1]['y']],
                        color=line_c, linewidth=pair_line_width, linestyle=pair_line_style, zorder=1,
                        alpha=pair_line_alpha)
                
    return ax
