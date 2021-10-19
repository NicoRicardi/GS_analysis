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

en = pd.read_csv("results_en_kcal.csv", index_col=0)
dens = pd.read_csv("results_densities.csv", index_col=0)

systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
syst_lbls = {"7HQ_2MeOH": "7HQ 2MeOH",
             "7HQ_formate": "7HQ formate",
             "Uracil_5H2O": "Uracil 5H2O",
             "XVI_2HCOOH": "XVI 2HCOOH"}

#r"$\int\frac{|\rho(\vec{r}|)}{2}d\vec{r}$"
en["E_MPk"] = en["E_FDET_MP"] - en["kernel_tot"]
s_ens, s_dens, en_ranges, en_centers = [], [], [], []
dens_rangeM, dens_cenM, dens_rangeP, dens_cenP = [], [], [], []
for n,system in enumerate(systems):
    s_en = en[en["system"] == system].set_index("calc")
    s_ens.append(s_en)
    s_den = dens[dens["system"] == system].set_index("calc")
    s_dens.append(s_den)
    all_ens = s_en["E_MPk"].tolist() + s_en["E_FDET_HF"].tolist() + [s_en["E_ref_HF"][0], s_en["E_ref_HF_CP"][0], s_en["E_ref_MP"][0], s_en["E_ref_MP_CP"][0]]
    en_ranges.append(max(all_ens) - min(all_ens))
    en_centers.append(0.5*(max(all_ens) + min(all_ens)))
    dens_rangeM.append(s_den["M_value"].max() - s_den["M_value"].min())
    dens_cenM.append(0.5*(s_den["M_value"].max() + s_den["M_value"].min()))
    dens_rangeP.append(s_den["densdiff_FDET_ref"].max() - s_den["densdiff_FDET_ref"].min())
    dens_cenP.append(0.5*(s_den["densdiff_FDET_ref"].max() + s_den["densdiff_FDET_ref"].min()))
    
en_range = max(en_ranges)*1.1   
den_rangeM = max(dens_rangeM)*1.1
den_rangeP = max(dens_rangeP)*1.1



colours = ["r", "g", "b", "c"]
axes_lbls = {"E_FDET_HF": r"$E^{FDET}_{HF} [Kcal/mol]$ ",
             "E_FDET_MP": r"$E^{FDET}_{MP,k}[Kcal/mol]$",
             "densdiff_FDET_ref": "P[a.u.]",
             "M_value": "M[a.u.]",
             "E_ref_HF": r"$E^{HF}_{}$",
             "E_ref_HF_CP": r"$E^{HF}_{CP}$",
             "E_ref_MP": r"$E^{MP}_{}$",
             "E_ref_MP_CP": r"$E^{MP}_{CP}$"}
lbls = ["none", "Mulliken", "ChelPG", "FT"]
fig_P, fig_HF, fig_PHF, fig_MP, fig_PMP = [plt.figure(figsize=(20, 10), dpi=150) for i in range(5)]
axs_P, axs_HF, axs_PHF, axs_MP, axs_PMP = [[] for i in range(5)]
constant = dict(linestyle="", alpha=0.75, markeredgecolor="k", markeredgewidth=2.5, markersize=15)

for n, system in enumerate(systems):
    axs_P.append(fig_P.add_subplot(221+n))#, label="P"))
    axs_HF.append(fig_HF.add_subplot(221+n))#, label="HF"))
    axs_PHF.append(fig_PHF.add_subplot(221+n))#, label="PHF"))
    axs_MP.append(fig_MP.add_subplot(221+n))#, label="MP"))
    axs_PMP.append(fig_PMP.add_subplot(221+n))#, label="MP"))
    for m in range(4):
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
        axs_MP[n].set_xlim([Mmin, Mmax])
        axs_MP[n].set_ylim([emin, emax])
        axs_MP[n].grid(b=True)
        axs_PHF[n].set_xlim([Pmin, Pmax])
        axs_PHF[n].set_ylim([emin, emax])
        axs_PHF[n].grid(b=True)
        axs_PMP[n].set_xlim([Pmin, Pmax])
        axs_PMP[n].set_ylim([emin, emax])
        axs_PMP[n].grid(b=True)
        if m == 0:
            axs_HF[n].axhline(y=s_ens[n]["E_ref_HF"][0], color="k", alpha=0.5, linestyle=":")
            axs_HF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF"][0], axes_lbls["E_ref_HF"], color="k", ha="left", va="top")
            axs_HF[n].axhline(y=s_ens[n]["E_ref_HF_CP"][0], color="k", alpha=0.5, linestyle="--")
            axs_HF[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_HF_CP"][0], axes_lbls["E_ref_HF_CP"], color="k", ha="left", va="bottom")
            axs_PHF[n].axhline(y=s_ens[n]["E_ref_HF"][0], color="k", alpha=0.5, linestyle=":")
            axs_PHF[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_HF"][0], axes_lbls["E_ref_HF"], color="k", ha="left", va="top")
            axs_PHF[n].axhline(y=s_ens[n]["E_ref_HF_CP"][0], color="k", alpha=0.5, linestyle="--")
            axs_PHF[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_HF_CP"][0], axes_lbls["E_ref_HF_CP"], color="k", ha="left", va="bottom")
            axs_MP[n].axhline(y=s_ens[n]["E_ref_MP"][0], color="k", alpha=0.5, linestyle=":")
            axs_MP[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_MP"][0], axes_lbls["E_ref_MP"], color="k", ha="left", va="top")
            axs_MP[n].axhline(y=s_ens[n]["E_ref_MP_CP"][0], color="k", alpha=0.5, linestyle="--")
            axs_MP[n].text(Mmax + 0.01*den_rangeM, s_ens[n]["E_ref_MP_CP"][0], axes_lbls["E_ref_MP_CP"], color="k", ha="left", va="bottom")
            axs_PMP[n].axhline(y=s_ens[n]["E_ref_MP"][0], color="k", alpha=0.5, linestyle=":")
            axs_PMP[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_MP"][0], axes_lbls["E_ref_MP"], color="k", ha="left", va="top")
            axs_PMP[n].axhline(y=s_ens[n]["E_ref_MP_CP"][0], color="k", alpha=0.5, linestyle="--")
            axs_PMP[n].text(Pmax + 0.01*den_rangeP, s_ens[n]["E_ref_MP_CP"][0], axes_lbls["E_ref_MP_CP"], color="k", ha="left", va="bottom")
        
        axs_P[n].plot(s_dens[n]["M_value"][m], s_dens[n]["densdiff_FDET_ref"][m],
           marker=(3+m, 0, 0), color=colours[m], label=lbls[m], **constant)
        axs_HF[n].plot(s_dens[n]["M_value"][m], s_ens[n]["E_FDET_HF"][m],
           marker=(3+m, 0, 0), color=colours[m], label=lbls[m], **constant)
        axs_PHF[n].plot(s_dens[n]["densdiff_FDET_ref"][m], s_ens[n]["E_FDET_HF"][m],
           marker=(3+m, 0, 0), color=colours[m], label=lbls[m], **constant)
        axs_MP[n].plot(s_dens[n]["M_value"][m], s_ens[n]["E_MPk"][m],
           marker=(3+m, 0, 0), color=colours[m], label=lbls[m], **constant)
        axs_PMP[n].plot(s_dens[n]["densdiff_FDET_ref"][m], s_ens[n]["E_MPk"][m],
           marker=(3+m, 0, 0), color=colours[m], label=lbls[m], **constant)
        axs_P[n].set_xlabel(axes_lbls["M_value"])
        axs_HF[n].set_xlabel(axes_lbls["M_value"])
        axs_PHF[n].set_xlabel(axes_lbls["densdiff_FDET_ref"])
        axs_MP[n].set_xlabel(axes_lbls["M_value"])
        axs_PMP[n].set_xlabel(axes_lbls["densdiff_FDET_ref"])
        axs_P[n].set_ylabel(axes_lbls["densdiff_FDET_ref"])
        axs_HF[n].set_ylabel(axes_lbls["E_FDET_HF"])
        axs_PHF[n].set_ylabel(axes_lbls["E_FDET_HF"])
        axs_MP[n].set_ylabel(axes_lbls["E_FDET_MP"])
        axs_PMP[n].set_ylabel(axes_lbls["E_FDET_MP"])
        axs_P[n].set_title(syst_lbls[system])
        axs_HF[n].set_title(syst_lbls[system])
        axs_PHF[n].set_title(syst_lbls[system])
        axs_MP[n].set_title(syst_lbls[system])
        axs_PMP[n].set_title(syst_lbls[system])
#        axs_P[n].
#        axs_HF[n].
#        axs_PHF[n].
#        axs_MP[n].
#        axs_PMP[n].
        if n == 0:
            axs_P[n].legend(loc="best")
            axs_HF[n].legend(loc="best")
            axs_PHF[n].legend(loc="best")
            axs_MP[n].legend(loc="best")
            axs_PMP[n].legend(loc="best")

fig_P.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
fig_HF.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
fig_PHF.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
fig_MP.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)
fig_PMP.subplots_adjust(top=0.96, bottom=0.065, left=0.055, right=0.97, wspace=0.18, hspace=0.25)

fig_P.savefig("M_vs_P.png")
fig_HF.savefig("M_vs_HF.png")
fig_PHF.savefig("P_vs_HF.png")
fig_MP.savefig("M_vs_MP.png")
fig_PMP.savefig("P_vs_MP.png")

###############
####PLOT ALL IN ONE PICTURE
###############
#for n, system in enumerate(systems):
#    axs.append(fig.add_subplot(221+n))
#    twins.append(axs[n].twinx())
#    for m in range(4):
#        axs[n].plot(s_dens[n]["M_value"][m], s_ens[n][props[0]][m], linestyle="",
#           marker=(3+m, 0, 0), color=colours[0], label=labels[0], alpha=0.5)
#        axs[n].plot(s_dens[n]["M_value"][m], s_ens[n]["E_MPk"][m], linestyle="",
#           marker=(3+m, 0, 0), color=colours[1], label=labels[1], alpha=0.5)
#        xmin = dens_cenx[n] - 0.5*den_rangex
#        xmin = max(0, xmin)
#        xmax = xmin + den_rangex
#        axs[n].set_xlim([xmin, xmax])
#        emin = en_centers[n] - 0.5*en_range
#        emax = emin + en_range
#        axs[n].set_ylim([emin, emax])
###        axs[n].legend()
#        twins[n].plot(s_dens[n]["M_value"][m], s_dens[n][props[2]][m], linestyle="",
#           marker=(3+m, 0, 0), color=colours[2], label=labels[2], alpha=0.5)
#        ymin = dens_ceny[n] - 0.5*den_rangey
#        ymax = ymin + den_rangey
#        twins[n].set_ylim([ymin - 0.05*den_rangey, ymax + 0.1*den_rangey])
##        twins[n].legend()
#fig.savefig("picture.png")