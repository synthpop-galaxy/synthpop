# Test script for Charon interpolator
# Authors: Jonas Klueter & Macy Huston

import sys
import os
sys.path.append(f'{os.path.abspath(os.path.dirname(__file__))}/..')
#print(sys.path)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from modules.evolution.mist import MIST
from modules.evolution.charon_interpolator import CharonInterpolator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

columns = ["EEP","UBVRIplus", "log10_isochrone_age_yr", "initial_mass", "[Fe/H]_init", "phase", 'star_mass']

mist = MIST(columns)

aa = mist.read_csv("pop5.cmd")
#aa = mist.read_csv("MIST_iso_63fea4a6b8c42.iso.cmd")

#isochrones = pd.concat(mist.isochrones.values)
grouped = mist.isochrones_grouped

mod = CharonInterpolator(isochrones=mist.isochrones)
mod2 = CharonInterpolator(isochrones=mist.isochrones, props_no_charon={"Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"})
age = np.linspace(7,10.3,100000)
"""
plt.semilogy(mod.tips.loc[-0.25].tip0,"r+")
plt.semilogy(mod.tips.loc[-0.25].tip1,'g+')
plt.semilogy(mod.tips.loc[-0.25].tip2,'b+')
plt.semilogy(age, mod.tips0_func(-0.25, age, grid=False), 'r')
plt.semilogy(age, mod.tips1_func(-0.25, age, grid=False), 'g')
plt.semilogy(age, mod.tips2_func(-0.25, age, grid=False), 'b')
plt.xlim([7,10])
plt.ylim([0.5,25])
"""

nn = 100_000
age = 10**(aa.log10_isochrone_age_yr[0]-9) * np.ones(nn)
met = aa['[Fe/H]_init'][0]* np.ones(nn)
min_mass1 = mod.tips1_func(met[0], np.log10(age[0])+9, grid=False)
min_mass2 = mod.tips2_func(met[0], np.log10(age[0])+9, grid=False)
imass1 =np.linspace(0.1,min_mass1, nn//6)
imass2 =np.linspace(min_mass1, min_mass2, nn//2)
imass3 = np.logspace(-10, np.log10(max(aa.initial_mass)-min_mass2+1), nn-nn//3*2) + min_mass2
     # np.linspace(min_mass2, max(aa.initial_mass), nn-nn//6-nn//2)
imass = np.hstack([imass1,imass2,imass3])

metlow,methigh,lower,upper, _, _ = mod.get_adjacent_gridpoints(met, np.log10(age)+9)
iso1 = mod.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"]).get_group((methigh[0], lower[0]))
iso2 = mod.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"]).get_group((methigh[0], upper[0]))
prop = "EEP"
output = mod.get_evolved_props(imass, met, age, {prop, "Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"})[0]
output2 = mist.get_evolved_props(imass, met, age,{prop, "Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"})[0]
"""
w1,w2,w3,w4 = mod.get_weights(met, np.log10(age)+9,metlow, methigh, lower,upper )
print(w1[0], w2[0], w3[0], w4[0])

mass1 = mod.get_modified_mass(imass,met, np.log10(age)+9,  methigh, lower)
mass2 = mod.get_modified_mass(imass,met, np.log10(age)+9,  methigh, upper)
pp = { prop, "Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"}
p_f1_charon, in_grid1 = mod.interp_props(mass1, methigh, lower, pp)
p_f2_charon, in_grid2 = mod.interp_props(mass2, methigh, upper, pp)

plt.plot(aa.initial_mass, aa[prop],'g')
iprop = list(pp).index(prop)
plt.plot(imass, p_f1_charon[:,iprop],'b')
plt.plot(mass1, p_f1_charon[:,iprop],'b--')
plt.plot(imass, p_f2_charon[:,iprop],'r')
plt.plot(mass2, p_f2_charon[:,iprop],'r--')


#plt.plot(imass, output[prop],'g')
#plt.plot(imass, output2[prop],'r')
plt.axvline(mod.endp_func(met[0], np.log10(age[0])+9, grid=False), color="k")
plt.axvline(mod.tips2_func(met[0], np.log10(age[0])+9, grid=False), color="k")
plt.axvline(mod.endp_func(met[0], lower[0], grid=False), color="k", linestyle=":")
plt.axvline(mod.tips2_func(met[0], lower[0], grid=False), color="k", linestyle=":")
plt.axvline(mod.endp_func(met[0], upper[0], grid=False), color="k", linestyle="--")
plt.axvline(mod.tips2_func(met[0], upper[0], grid=False), color="k",linestyle="--")

plt.figure()
bb1 = mod.grouped_iso.get_group((methigh[0], lower[0]))
bb2 = mod.grouped_iso.get_group((methigh[0], upper[0]))

plt.plot(bb1.EEP, bb1[prop],'r')
plt.plot(bb2.EEP, bb2[prop],'g')
plt.plot(aa.EEP, aa[prop],'b')



ig = list(pp).index("Gaia_G_EDR3")
irp = list(pp).index("Gaia_RP_EDR3")
ibp = list(pp).index("Gaia_BP_EDR3")
plt.figure()


plt.plot(bb2.Gaia_BP_EDR3-bb2.Gaia_RP_EDR3, bb2.Gaia_G_EDR3,'k', lw=0.5)

plt.plot(aa.Gaia_BP_EDR3-aa.Gaia_RP_EDR3, aa.Gaia_G_EDR3,'g', lw=8, alpha=0.2)
plt.plot(aa.Gaia_BP_EDR3-aa.Gaia_RP_EDR3, aa.Gaia_G_EDR3,'g', lw=8, alpha=0.2)
plt.plot(aa.Gaia_BP_EDR3-aa.Gaia_RP_EDR3, aa.Gaia_G_EDR3,'g', lw=8, alpha=0.2)
plt.plot(output["Gaia_BP_EDR3"]-output["Gaia_RP_EDR3"], output["Gaia_G_EDR3"],'g.',ms=1)
plt.plot(output2["Gaia_BP_EDR3"]-output2["Gaia_RP_EDR3"], output2["Gaia_G_EDR3"],'r.', ms=1)
plt.plot(p_f1_charon[:,ibp]-p_f1_charon[:,irp], p_f1_charon[:,ig],'c', lw=0.5)
plt.plot(p_f2_charon[:,ibp]-p_f2_charon[:,irp], p_f2_charon[:,ig],'b', lw=0.5)
print(len(aa))
print(max(output["Gaia_BP_EDR3"]))
print(len(output["Gaia_BP_EDR3"]))
plt.gca().invert_yaxis()


"""
def figcmd(ax, G, BP, RP, G1, BP1, RP1):
    ax.plot( BP - RP, G, color="b", label="CharonInterpolator", zorder=4,)
    ax.plot( BP1 - RP1, G1, color="r", label="LagrangeInterpolator", zorder=3,)

    ax.plot(aa.Gaia_BP_EDR3 - aa.Gaia_RP_EDR3, aa.Gaia_G_EDR3, 'g', lw=8, alpha=0.2, label= "MIST Web Interpolator")

    ax.set_xlabel("BP - RP")
    ax.set_ylabel("G")
    ax.legend(loc=4)
    ax.invert_yaxis()
    ax.set_xlim(-1.2,4.9)
    #ax.set_ylim(12.5, -5.0)

def figimass(ax,mod, imass,G1, G2, iso1, iso2):

    ax.plot(imass, G1, "b", zorder=4, label="CharonInterpolator")
    ax.plot(imass, G2, "r", zorder=4, lw=2, label="LagrangeInterpolator")
    ax.plot(aa.initial_mass, aa.Gaia_G_EDR3, 'g', lw=8, alpha=0.2, label = "MIST Web Interpolator")

    ax.plot(iso1.initial_mass, iso1.Gaia_G_EDR3, 'k--', label="Adjacent grid points",zorder=5)
    ax.plot(iso2.initial_mass, iso2.Gaia_G_EDR3, 'k--',zorder=5)

    ax.set_xlabel(r"m$_{\rm init}$ (M$_\odot$)")
    #ax.set_ylabel("G")
    ax.invert_yaxis()
    ax.legend(loc=4)
    ax.set_xlim(0.952, 1.025)
    
def figimass_zoom(ax,mod, imass,G, iso1, iso2):
    tip = mod.tips2_func(aa["[Fe/H]_init"], aa["log10_isochrone_age_yr"], grid=False)[0]
    offset = -tip
    ax.plot(imass+offset, G, "b", zorder=4, label="CharonInterpolator")
    ax.plot(aa.initial_mass+offset, aa.Gaia_G_EDR3, 'g',lw=8, alpha=0.2, label = "MSIT Web Interpolator")

    modmass1 = mod.get_modified_mass(iso1.initial_mass,
            iso1["[Fe/H]_init"], iso1["log10_isochrone_age_yr"],
            aa["[Fe/H]_init"].iloc[0]*np.ones(len(iso1.initial_mass)), aa["log10_isochrone_age_yr"].iloc[0]*np.ones(len(iso1.initial_mass)),)
    ax.plot(modmass1+offset, iso1.Gaia_G_EDR3, 'k:', label="CharonInterpolator\n shifted grid points")

    modmass2 = mod.get_modified_mass(iso2.initial_mass,
            iso2["[Fe/H]_init"], iso2["log10_isochrone_age_yr"],
            aa["[Fe/H]_init"].iloc[0]*np.ones(len(iso2.initial_mass)), aa["log10_isochrone_age_yr"].iloc[0]*np.ones(len(iso2.initial_mass)),)
    #tip_iso = aa.initial_mass.loc[aa.EEP == 1409].iloc[0]

    ax.plot(modmass2+offset, iso2.Gaia_G_EDR3, 'k:')
    ax.axvline(tip+offset,color="k", ls="--", label = "Est. end of AGB")
    #ax.axvline(tip_iso,color="r", label = "end of the AGB from isochrones ")
    ax.invert_yaxis()
    ax.set_xlabel(r"m$_{\rm init}$-m$_{\rm ABG}$ (10$^{-6}$ M$_\odot$)")
    #ax.set_ylabel("G")
    ax.legend(loc=4)
    ax.set_xlim(0.9938135+offset, 0.9938163+offset)
    #ax.set_xscale('log')
    #ax.figure.canvas.draw()
    #offset = ax.xaxis.get_major_formatter().get_offset()
    #print(offset)
    ax.xaxis.offsetText.set_visible(False)
    #ax.xaxis.set_label_text("initial_mass" + " " + offset)

fig, axes = plt.subplots(1,3, figsize=(12,6), sharey = True)
print(axes)
figcmd(axes[0], output["Gaia_G_EDR3"], output["Gaia_BP_EDR3"], output["Gaia_RP_EDR3"],  \
       output2["Gaia_G_EDR3"],output2["Gaia_BP_EDR3"], output2["Gaia_RP_EDR3"])
figimass(axes[1],mod, imass, output["Gaia_G_EDR3"], output2["Gaia_G_EDR3"], iso1, iso2)
figimass_zoom(axes[2], mod, imass, output["Gaia_G_EDR3"], iso1, iso2)
fig.subplots_adjust(wspace=0,bottom=0.2, left=0.1)

plt.savefig('validation_figures/charon_fig_test.pdf')

