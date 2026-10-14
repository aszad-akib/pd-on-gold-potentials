#!/usr/bin/env python
"""Pd adsorption on Au(111), Au(100), Au(110): every high-symmetry site, one potential.
Built from Pt_on_gold/Pt_on_Au-slab_EMT.py in miroi/open-collection, with these changes:
  - each potential uses its OWN lattice constant (equation-of-state fit)
  - the clean slab is relaxed with the same constraints as the slab + Pd
  - bottom half of the slab frozen, top half and Pd relaxed
  - Pd is held laterally on its site (FixCartesian), so each energy belongs to that site
usage:  python pd_production.py lj|emt|eam|mace|mpa|chgnet     ->  results_<pot>.csv
"""
import sys, csv, time, warnings
import numpy as np
from ase import Atoms
from ase.build import bulk, fcc111, fcc100, fcc110, add_adsorbate
from ase.constraints import FixAtoms, FixCartesian
from ase.eos import EquationOfState
from ase.optimize import BFGS
from potentials import get_calc
warnings.filterwarnings("ignore")

pot = sys.argv[1] if len(sys.argv) > 1 else "emt"
calc = get_calc(pot)
t0 = time.time()

# 1. this potential's own lattice constant of gold
vols, ens = [], []
for a in np.linspace(3.90, 4.30, 9):
    au = bulk("Au", "fcc", a=a, cubic=True); au.calc = calc
    vols.append(au.get_volume()); ens.append(au.get_potential_energy())
a0 = EquationOfState(vols, ens).fit()[0] ** (1 / 3)

# 2. isolated Pd atom (same potential)
pd = Atoms("Pd", cell=[18, 18, 18], pbc=True); pd.center(); pd.calc = calc
e_pd = pd.get_potential_energy()
print(f"[{pot}] a0 = {a0:.3f} A (exp 4.078)   E(Pd atom) = {e_pd:.3f} eV", flush=True)

facets = {
    "Au(111)": (fcc111, 4, ["ontop", "bridge", "fcc", "hcp"]),
    "Au(100)": (fcc100, 4, ["ontop", "bridge", "hollow"]),
    "Au(110)": (fcc110, 6, ["ontop", "shortbridge", "longbridge", "hollow"]),
}
rows = []
for facet, (build, layers, sites) in facets.items():
    slab = build("Au", size=(3, 3, layers), a=a0, vacuum=12.0)
    slab.set_constraint(FixAtoms(mask=[t > layers // 2 for t in slab.get_tags()]))
    slab.calc = calc
    BFGS(slab, logfile=None).run(fmax=0.02, steps=500)
    e_slab = slab.get_potential_energy()
    ztop = slab.positions[:, 2].max()
    for site in sites:
        s = slab.copy()
        add_adsorbate(s, "Pd", 2.0, site)
        s.set_constraint(s.constraints + [FixCartesian(len(s) - 1, mask=(1, 1, 0))])
        s.calc = calc
        BFGS(s, logfile=None).run(fmax=0.02, steps=800)
        e = s.get_potential_energy() - e_slab - e_pd
        d = s.get_distances(len(s) - 1, range(len(s) - 1), mic=True)
        nb = int((d < d.min() + 0.3).sum())
        h = s.positions[-1, 2] - ztop
        rows.append(dict(potential=pot, facet=facet, site=site, n_Au_neighbours=nb,
                         E_ads_eV=round(e, 4), E_ads_kJmol=round(e * 96.485, 1),
                         height_A=round(h, 2), d_PdAu_A=round(d.min(), 3)))
        print(f"{facet:8s} {site:12s} {e:8.3f} eV {e * 96.485:8.1f} kJ/mol  h={h:5.2f}  d(Pd-Au)={d.min():.2f}", flush=True)
        s.write(f"Pd_{facet[3:6]}_{site}_{pot}.xyz")

with open(f"results_{pot}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"saved results_{pot}.csv   ({time.time() - t0:.0f} s)")
