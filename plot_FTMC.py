#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:52:30 2021

@author: nico
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 08:31:39 2021

@author: nico
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', size=14)
plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=14)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('figure', titlesize=14)
plt.rc('legend', fontsize=14)


def get_iter(string):
    num = [i for i in string if i.isnumeric()]
    if num:
        return int("".join(num))
    else:
        return -1

raw_en = pd.read_csv("FT_en_kcal.csv", index_col=0)
raw_dens = pd.read_csv("FT_dens.csv", index_col=0)
raw_MCen = pd.read_csv("FTMC_en_kcal.csv", index_col=0)
raw_MCdens = pd.read_csv("FTMC_dens.csv", index_col=0)

systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
syst_lbls = {"7HQ_2MeOH": "7HQ 2MeOH",
             "7HQ_formate": "7HQ formate",
             "Uracil_5H2O": "Uracil 5H2O",
             "XVI_2HCOOH": "XVI 2HCOOH"}

#raw_en["E_MPk"] = raw_en["E_FDET_MP"] - raw_en["kernel_tot"]

colours = ["r", "g", "b", "c"]
axes_lbls = {"E_FDET_HF": r"$E^{FDET}_{HF} [Kcal/mol]$ ",
             "E_FDET_MP": r"$E^{FDET}_{MP,k}[Kcal/mol]$",
             "densdiff_FDET_ref": "P[a.u.]",
             "M_value": "M[a.u.]",
             "E_ref_HF": r"$E^{HF}_{}$",
             "E_ref_HF_CP": r"$E^{HF}_{CP}$",
             "E_ref_MP": r"$E^{MP}_{}$",
             "E_ref_MP_CP": r"$E^{MP}_{CP}$"}


constant = dict(linestyle="", alpha=0.35, markeredgecolor="k", markeredgewidth=2.5, markersize=2.5)

for exp in ["ME", "SE"]:
#for exp in ["ME"]:
    en = raw_en[raw_en["calc"].apply(lambda x: exp in x)]
    MCen = raw_MCen[raw_MCen["calc"].apply(lambda x: exp in x)]
    dens = raw_dens[raw_dens["calc"].apply(lambda x: exp in x)]
    MCdens = raw_MCdens[raw_MCdens["calc"].apply(lambda x: exp in x)]
    s_ens, s_MCens, s_dens, s_MCdens, en_ranges, en_centers = [], [], [], [], [], []
    dens_rangeM, dens_cenM, dens_rangeP, dens_cenP = [], [], [], []
    for n,system in enumerate(systems):
        s_en = en[en["system"] == system]
        s_en.loc[:,"iter"] = s_en["calc"].apply(get_iter)
        s_en.loc[s_en["iter"] == -1,"iter"] = s_en["iter"].max() + 2
        s_en = s_en.set_index("iter").sort_index()#[1:]
        s_ens.append(s_en)
        s_den = dens[dens["system"] == system]
        s_den.loc[:,"iter"] = s_den["calc"].apply(get_iter)
        s_den.loc[s_den["iter"] == -1,"iter"] = s_den["iter"].max() + 2
        s_den = s_den.set_index("iter").sort_index()#[1:]
        s_dens.append(s_den)
        s_MCen = MCen[MCen["system"] == system]
        s_MCen.loc[:,"iter"] = s_MCen["calc"].apply(get_iter)
        s_MCen.loc[s_MCen["iter"] == -1,"iter"] = s_MCen["iter"].max() + 2
        s_MCen = s_MCen.set_index("iter").sort_index()#[1:]
        s_MCens.append(s_MCen)
        s_MCden = MCdens[MCdens["system"] == system]
        s_MCden.loc[:,"iter"] = s_MCden["calc"].apply(get_iter)
        s_MCden.loc[s_MCden["iter"] == -1,"iter"] = s_MCden["iter"].max() + 2
        s_MCden = s_MCden.set_index("iter").sort_index()#[1:]
        s_MCdens.append(s_MCden)
        all_ens = s_en["E_FDET_HF"].tolist() + [s_en["E_ref_HF"][2], s_en["E_ref_HF_CP"][2]] + s_MCen["E_FDET_HF"].tolist() + [s_MCen["E_ref_HF"][2], s_MCen["E_ref_HF_CP"][2]]
        en_ranges.append(max(all_ens) - min(all_ens))
        en_centers.append(0.5*(max(all_ens) + min(all_ens)))
        dens_rangeM.append(s_den["M_value"].max() - s_den["M_value"].min())
        dens_cenM.append(0.5*(s_den["M_value"].max() + s_den["M_value"].min()))
        dens_rangeP.append(max(s_den["densdiff_FDET_ref"].max(), s_MCden["densdiff_FDET_ref"].max()) - min(s_den["densdiff_FDET_ref"].min(), s_MCden["densdiff_FDET_ref"].min()))
        dens_cenP.append(0.5*(max(s_den["densdiff_FDET_ref"].max(), s_MCden["densdiff_FDET_ref"].max()) + min(s_den["densdiff_FDET_ref"].min(), s_MCden["densdiff_FDET_ref"].min())))
    
    en_range = max(en_ranges)*1.1   
    den_rangeM = max(dens_rangeM)*1.1
    den_rangeP = max(dens_rangeP)*1.1
    fig_P, fig_HF, fig_PHF = [plt.figure(figsize=(20, 10), dpi=150) for i in range(3)]
    axs_P, axs_HF, axs_PHF = [[] for i in range(3)]
    fig_NM, fig_NP, fig_NHF = [plt.figure(figsize=(20, 10), dpi=150) for i in range(3)]
    axs_NM, axs_NP, axs_NHF = [[] for i in range(3)]
    for n, system in enumerate(systems):
        axs_P.append(fig_P.add_subplot(221+n))#, label="P"))
        axs_HF.append(fig_HF.add_subplot(221+n))#, label="HF"))
        axs_PHF.append(fig_PHF.add_subplot(221+n))#, label="PHF"))
        axs_NM.append(fig_NM.add_subplot(221+n))#, label="NM"))
        axs_NP.append(fig_NP.add_subplot(221+n))#, label="NP"))
        axs_NHF.append(fig_NHF.add_subplot(221+n))#, label="NHF"))
        Mmin = dens_cenM[n] - 0.5*den_rangeM
        Mmin = max(0, Mmin)
        Mmax = Mmin + den_rangeM
        emin = en_centers[n] - 0.5*en_range
        emax = emin + en_range
        Pmin = dens_cenP[n] - 0.5*den_rangeP
        Pmin = max(Pmin, 0) 
        Pmax = Pmin + den_rangeP
        axs_P[n].set_xlim([Mmin, Mmax])
        axs_P[n].set_ylim([Pmin, Pmax])
        axs_P[n].grid(b=True)
        axs_HF[n].set_xlim([Mmin, Mmax])
        axs_HF[n].set_ylim([emin, emax])
        axs_HF[n].grid(b=True)
        axs_PHF[n].set_xlim([Pmin, Pmax])
        axs_PHF[n].set_ylim([emin, emax])
        axs_PHF[n].grid(b=True)
        axs_NM[n].set_ylim([Mmin, Mmax])
        axs_NM[n].grid(b=True)
        axs_NP[n].set_ylim([Pmin, Pmax])
        axs_NP[n].grid(b=True)
        axs_NHF[n].set_ylim([emin, emax])
        axs_NHF[n].grid(b=True)
        axs_HF[n].axhline(y=s_ens[n]["E_ref_HF"][2], color="k", alpha=0.5, linestyle=":")
        axs_HF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF"][2], axes_lbls["E_ref_HF"], color="k", ha="left", va="top")
        axs_HF[n].axhline(y=s_ens[n]["E_ref_HF_CP"][2], color="k", alpha=0.5, linestyle="--")
        axs_HF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF_CP"][2], axes_lbls["E_ref_HF_CP"], color="k", ha="left", va="bottom")
        axs_NHF[n].axhline(y=s_ens[n]["E_ref_HF"][2], color="k", alpha=0.5, linestyle=":")
        axs_NHF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF"][2], axes_lbls["E_ref_HF"], color="k", ha="left", va="top")
        axs_NHF[n].axhline(y=s_ens[n]["E_ref_HF_CP"][2], color="k", alpha=0.5, linestyle="--")
        axs_NHF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF_CP"][2], axes_lbls["E_ref_HF_CP"], color="k", ha="left", va="bottom")
        axs_PHF[n].axhline(y=s_ens[n]["E_ref_HF"][2], color="k", alpha=0.5, linestyle=":")
        axs_PHF[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_HF"][2], axes_lbls["E_ref_HF"], color="k", ha="left", va="top")
        axs_PHF[n].axhline(y=s_ens[n]["E_ref_HF_CP"][2], color="k", alpha=0.5, linestyle="--")
        axs_PHF[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_HF_CP"][2], axes_lbls["E_ref_HF_CP"], color="k", ha="left", va="bottom")
        axs_P[n].plot(s_dens[n]["M_value"], s_dens[n]["densdiff_FDET_ref"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_P[n].plot(s_MCdens[n]["M_value"], s_dens[n]["densdiff_FDET_ref"],
           marker="o", color=colours[2], label="_nolegend_", **constant)
        axs_HF[n].plot(s_dens[n]["M_value"], s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_HF[n].plot(s_MCdens[n]["M_value"], s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[2], label="_nolegend_", **constant)
        axs_PHF[n].plot(s_dens[n]["densdiff_FDET_ref"], s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_PHF[n].plot(s_MCdens[n]["densdiff_FDET_ref"], s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[2], label="_nolegend_", **constant)
        axs_NM[n].plot(s_dens[n].index, s_dens[n]["M_value"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_NP[n].plot(s_dens[n].index, s_dens[n]["densdiff_FDET_ref"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_NP[n].plot(s_MCdens[n].index, s_dens[n]["densdiff_FDET_ref"],
           marker="o", color=colours[2], label="_nolegend_", **constant)
        axs_NHF[n].plot(s_dens[n].index, s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[0], label="_nolegend_", **constant)
        axs_NHF[n].plot(s_MCdens[n].index, s_ens[n]["E_FDET_HF"],
           marker="o", color=colours[2], label="_nolegend_", **constant)
        axs_P[n].set_xlabel(axes_lbls["M_value"])
        axs_HF[n].set_xlabel(axes_lbls["M_value"])
        axs_PHF[n].set_xlabel(axes_lbls["densdiff_FDET_ref"])
        axs_P[n].set_ylabel(axes_lbls["densdiff_FDET_ref"])
        axs_HF[n].set_ylabel(axes_lbls["E_FDET_HF"])
        axs_PHF[n].set_ylabel(axes_lbls["E_FDET_HF"])
        axs_NM[n].set_ylabel(axes_lbls["M_value"])
        axs_NP[n].set_ylabel(axes_lbls["densdiff_FDET_ref"])
        axs_NHF[n].set_ylabel(axes_lbls["E_FDET_HF"])
        axs_P[n].set_title(syst_lbls[system])
        axs_HF[n].set_title(syst_lbls[system])
        axs_PHF[n].set_title(syst_lbls[system])
        axs_NM[n].set_title(syst_lbls[system])
        axs_NP[n].set_title(syst_lbls[system])
        axs_NHF[n].set_title(syst_lbls[system])
#        axs_NM[n].set_yscale("symlog")
#        axs_NP[n].set_yscale("symlog")
#        axs_NHF[n].set_yscale("symlog")
#        axs_NMP[n].set_yscale("symlog")
    fig_P.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    fig_HF.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    fig_PHF.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    fig_NM.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    fig_NP.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    fig_NHF.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
    
    fig_P.savefig("{}_cycles_M_vs_P.png".format(exp))
    fig_HF.savefig("{}_cycles_M_vs_HF.png".format(exp))
    fig_PHF.savefig("{}_cycles_P_vs_HF.png".format(exp))
    fig_NM.savefig("{}_cycles_N_vs_M.png".format(exp))
    fig_NP.savefig("{}_cycles_N_vs_P.png".format(exp))
    fig_NHF.savefig("{}_cycles_N_vs_HF.png".format(exp))


