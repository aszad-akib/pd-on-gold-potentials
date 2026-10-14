#!/usr/bin/env python
"""Figures for Pd on gold:  figures/Eads_ladder.png  and  figures/Pd_alloying.png
reads results/results_<pot>.csv and results/alloy_<pot>.csv"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POTS = [("lj", "Lennard-Jones"), ("emt", "EMT"), ("eam", "EAM (Zhou 2004)"),
        ("mace", "MACE-MP-0"), ("mpa", "MACE-MPA-0"), ("chgnet", "CHGNet")]
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
# GGA-DFT, Gotsis et al., Chem. Phys. Lett. 442, 105 (2007)
DFT = {("Au(111)", "ontop"): -1.688, ("Au(111)", "fcc"): -2.564, ("Au(111)", "hcp"): -2.510,
       ("Au(100)", "ontop"): -1.943, ("Au(100)", "bridge"): -2.223, ("Au(100)", "hollow"): -2.995}
os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 10})

R = {p: list(csv.DictReader(open(f"results/results_{p}.csv"))) for p, _ in POTS}
sites = [(r["facet"], r["site"]) for r in R["emt"]]
E = {p: [float(r["E_ads_eV"]) for r in R[p]] for p, _ in POTS}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2.4, 1]})
w = 0.13
x = np.arange(len(sites))
for k, ((p, name), c) in enumerate(zip(POTS, COL)):
    ax1.bar(x + (k - 2.5) * w, E[p], w * 0.9, color=c, label=name)
dx = [i for i, s in enumerate(sites) if s in DFT]
ax1.scatter(dx, [DFT[sites[i]] for i in dx], marker="D", s=40, color="k", zorder=5, label="DFT (Gotsis 2007)")
ax1.set_xticks(x); ax1.set_xticklabels([f"{f[3:6]}\n{s}" for f, s in sites], fontsize=8)
ax1.set_ylabel("E_ads (eV)"); ax1.axhline(0, color="0.6", lw=0.8)
ax1.set_title("Pd adsorption on gold: 6 potentials, 3 facets, 11 sites")
ax1.legend(fontsize=8, ncol=4, loc="lower left", frameon=False)
ax1.grid(axis="y", color="0.9"); ax1.set_axisbelow(True)

ref = np.array([DFT[s] for s in sites if s in DFT])
for (p, name), c in zip(POTS, COL):
    y = np.array([E[p][i] for i in dx])
    mae = np.abs(y - ref).mean()
    ax2.scatter(ref, y, color=c, s=50, label=f"{name}  MAE {mae:.2f} eV", edgecolor="white", lw=1, zorder=3)
lim = [-4.0, 0.2]
ax2.plot(lim, lim, color="0.5", lw=1, ls="--"); ax2.set_xlim(lim); ax2.set_ylim(lim)
ax2.set_xlabel("DFT E_ads (eV)"); ax2.set_ylabel("potential E_ads (eV)")
ax2.set_title("Against DFT (6 sites on 111 and 100)")
ax2.legend(fontsize=7.5, loc="lower right", frameon=False)
plt.tight_layout(); plt.savefig("figures/Eads_ladder.png", dpi=150)
print("saved figures/Eads_ladder.png")

A = {p: next(csv.DictReader(open(f"results/alloy_{p}.csv"))) for p, _ in POTS}
fig, ax = plt.subplots(figsize=(8, 4.2))
lab = ["Pd swaps into top layer\n(exchange with an Au)", "Pd in layer 2\nvs layer 1", "Pd in layer 3\nvs layer 1"]
keys = ["dE_exchange_eV", "dE_L2_minus_L1_eV", "dE_L3_minus_L1_eV"]
x = np.arange(3)
for k, ((p, name), c) in enumerate(zip(POTS, COL)):
    ax.bar(x + (k - 2.5) * w, [float(A[p][q]) for q in keys], w * 0.9, color=c, label=name)
ax.axhline(0, color="0.4", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel("ΔE (eV)   (< 0: Pd goes into the gold)")
ax.set_title("Does Pd stay on the Au(111) surface?")
ax.legend(fontsize=8, ncol=2, frameon=False, loc="lower right")
ax.grid(axis="y", color="0.9"); ax.set_axisbelow(True)
plt.tight_layout(); plt.savefig("figures/Pd_alloying.png", dpi=150)
print("saved figures/Pd_alloying.png")

print(f"\n{'potential':16s} {'111 fcc':>8s} {'100 hol':>8s} {'110 hol':>8s} {'MAE vs DFT':>10s}")
for p, name in POTS:
    y = np.array([E[p][i] for i in dx])
    get = lambda f, s: E[p][sites.index((f, s))]
    print(f"{name:16s} {get('Au(111)','fcc'):8.2f} {get('Au(100)','hollow'):8.2f} "
          f"{get('Au(110)','hollow'):8.2f} {np.abs(y - ref).mean():10.2f}")
print(f"{'DFT (Gotsis)':16s} {-2.564:8.2f} {-2.995:8.2f} {'n/a':>8s}")
