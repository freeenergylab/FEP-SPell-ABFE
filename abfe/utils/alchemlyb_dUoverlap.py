#!/usr/bin/env python
# coding: utf-8
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This alchemlyb dUoverlap module is used to check and plot the convergence and
phase space overlap for free energy calculations.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (8, 6)
import pandas as pd
import seaborn as sns
from collections import defaultdict

from alchemlyb import concat
from alchemlyb.estimators import MBAR, BAR, TI
from alchemlyb.parsing import amber
from alchemlyb.parsing.amber import extract_u_nk, extract_dHdl
from alchemlyb.visualisation import plot_convergence
from alchemlyb.visualisation.dF_state import plot_dF_state
from alchemlyb.visualisation.ti_dhdl import plot_ti_dhdl
from alchemlyb.visualisation.mbar_matrix import plot_mbar_overlap_matrix

def beta(T):
    """Convert to beta in units of kT (beta in 1./(kcal/mol)).

    Args:
        T (float):
            temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT)

    Returns:
        beta: 1./(amber.k_b*T)
    """

    return 1./(amber.k_b*T)

def get_delta_f(estimator):
    """Get dF/ddF from estimator.

    Args:
        estimator (object): TI/BAR/MBAR fit()

    Returns:
        dF/ddF
    """
    ee = 0.0
    for i in range(len(estimator.d_delta_f_) - 1):
        ee += estimator.d_delta_f_.values[i][i+1]**2

    return estimator.delta_f_.iloc[0, -1], ee**0.5

def asserter(method, T):
    """The asserter is used for TI/BAR/MBAR fit().

    Args:
        methods (list):
            it contains all analysis methods needed to be done
        T (float):
            temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT)
    """
    assert method.delta_f_.attrs['temperature'] == T
    assert method.delta_f_.attrs['energy_unit'] == 'kT'
    assert method.d_delta_f_.attrs['temperature'] == T
    assert method.d_delta_f_.attrs['energy_unit'] == 'kT'
    if method == 'TI':
        assert method.dhdl.attrs['temperature'] == T
        assert method.dhdl.attrs['energy_unit'] == 'kT'

def plot_distribution_deltaU(u_nk, T, units):
    """Plot deltaU overlap distribution between neighbour windows.

    Args:
        u_nk (DataFrame):
            reduced potential for each alchemical state (k) for each frame (n).
        T (float):
            temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT)
        units (str):
            amber energy in kcal/mol
    """
    # sort by state so that rows from same state are in contiguous blocks
    u_nk = u_nk.sort_index(level=u_nk.index.names[1:])
    # get a list of the lambda states
    states_ = u_nk.columns.values.tolist()
    # group u_nk by lambda states
    groups = u_nk.groupby(level=u_nk.index.names[1:])
    N_k = [(len(groups.get_group(i)) if i in groups.groups else 0)
           for i in u_nk.columns]

    fig, axs = plt.subplots(len(N_k)-1, 1, figsize=(len(N_k), len(N_k)-1))
    fig.add_subplot(111, frameon=False)
    fig.tight_layout()
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)

    # Now get free energy differences and their uncertainties for each step
    for k in range(0, len(N_k) - 1, 1):
        # get us from lambda step k
        uk = groups.get_group(states_[k])
        # get w_F
        w_f = uk.iloc[:, k+1] - uk.iloc[:, k]
        # get us from lambda step k+1
        uk1 = groups.get_group(states_[k+1])
        # get w_R
        w_r = uk1.iloc[:, k] - uk1.iloc[:, k+1]
        sns.kdeplot(data=-1.0*w_r/beta(T), shade=True, legend=True,
                    label="%s-%s" % (str(k), str(k+1)), ax=axs[k],
                    color='cornflowerblue')
        sns.kdeplot(data=w_f/beta(T), shade=True, legend=True,
                    label="%s-%s" % (str(k+1), str(k)), ax=axs[k],
                    color='cornflowerblue')
        axs[k].legend(loc='best', frameon=False, fontsize=5)
        axs[k].set_xlabel(r'', fontsize=16)
        axs[k].set_ylabel(r'', fontsize=16)
    plt.title('Distribution of $\Delta $U=U(target)-U(reference) Overlap')
    plt.xlabel(rf'$\Delta U$ ({units})', fontsize=16)
    plt.ylabel(r'Density', fontsize=16)
    plt.savefig('distribution_deltaU.svg')

def alchemlyb_analysis(ti_out_files, methods, T, units):
    """Perform alchemical analysis and plots using alchemlyb.

    Args:
        ti_out_files (list):
            it contains all ti.out files
        methods (list):
            it contains all analysis methods needed to be done
        T (float):
            temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT)
        units (str):
            amber energy in kcal/mol

    Returns:
        df_all (DataFrame):
            it contains the dF/ddF for all analysis methods
    """
    data = defaultdict(list)
    data['dF'] = [None]*len(methods)
    data['ddF'] = [None]*len(methods)
    df_all = pd.DataFrame(data, index=methods)

    if 'MBAR' in methods or 'BAR' in methods:
        u_nk = concat([extract_u_nk(file, T=T) for file in ti_out_files])
        plot_distribution_deltaU(u_nk, T=T, units=units)
        if 'BAR' in methods:
            bar = BAR().fit(u_nk)
            asserter(bar, T)
            bar_delta_f, bar_d_delta_f = get_delta_f(bar)
            df_all['dF']['BAR'], df_all['ddF']['BAR'] = \
                bar_delta_f/beta(T), bar_d_delta_f/beta(T)
        if 'MBAR' in methods:
            mbar = MBAR().fit(u_nk)
            asserter(mbar, T)
            mbar_delta_f, mbar_d_delta_f = \
                mbar.delta_f_.iloc[0, -1], mbar.d_delta_f_.iloc[0, -1]
            df_all['dF']['MBAR'], df_all['ddF']['MBAR'] = \
                mbar_delta_f/beta(T), mbar_d_delta_f/beta(T)
            ax1 = plot_mbar_overlap_matrix(mbar.overlap_matrix)
            ax1.figure.savefig('mbar_overlap_matrix.svg')

            data_list = [extract_u_nk(file, T=T) for file in ti_out_files]
            forward = []
            forward_error = []
            backward = []
            backward_error = []
            num_points = 10
            for i in range(1, num_points+1):
                # Do the forward
                slic = int(len(data_list[0])/num_points*i)
                u_nk_forw = concat([data[:slic] for data in data_list])
                estimate = MBAR().fit(u_nk_forw)
                forward.append(estimate.delta_f_.iloc[0, -1])
                forward_error.append(estimate.d_delta_f_.iloc[0, -1])
                # Do the backward
                u_nk_back = concat([data[-slic:] for data in data_list])
                estimate = MBAR().fit(u_nk_back)
                backward.append(estimate.delta_f_.iloc[0, -1])
                backward_error.append(estimate.d_delta_f_.iloc[0, -1])
            forward = [value/beta(T) for value in forward]
            forward_error = [value/beta(T) for value in forward_error]
            backward = [value/beta(T) for value in backward]
            backward_error = [value/beta(T) for value in backward_error]
            try: 
                #v2.3.1
                data = {'Forward': forward, 'Forward_Error': forward_error,
                        'Backward': backward, 'Backward_Error': backward_error}
                df = pd.DataFrame(data)
                df.attrs['temperature'] = T
                df.attrs['energy_unit'] = units
                ax2 = plot_convergence(df, units=units)
            except TypeError:
                #v0.5.0
                ax2 = plot_convergence(forward, forward_error, backward, backward_error, units=units)
            ax2.figure.savefig('convergence.svg')

    if 'TI' in methods:
        dhdl = concat([extract_dHdl(file, T=T) for file in ti_out_files])
        ti = TI().fit(dhdl)
        asserter(ti, T)
        ti_delta_f = ti.delta_f_.iloc[0, -1]
        ti_d_delta_f = ti.d_delta_f_.iloc[0, -1]
        df_all['dF']['TI'], df_all['ddF']['TI'] = \
            ti_delta_f/beta(T), ti_d_delta_f/beta(T)
        dhdl_list = ti.separate_dhdl()
        for dhdl in dhdl_list:
            dhdl.attrs['temperature'] = T
            dhdl.attrs['energy_unit'] = 'kT'
        ax3 = plot_ti_dhdl(ti, units=units)
        ax3.figure.savefig('ti_dhdl.svg')

    estimators = []
    if 'TI' in methods:
        estimators.append(ti)
    if 'BAR' in methods:
        estimators.append(bar)
    if 'MBAR' in methods:
        estimators.append(mbar)
    # 'landscape' also can be treat as orientation
    ax4 = plot_dF_state(estimators, orientation='portrait', units=units)
    ax4.figure.savefig('dF_state.svg')

    return df_all
