#!/usr/bin/env python
"""Does a Pd atom stay on top of Au(111), or sink into the gold?  (3x3 cell, 5 layers)
  exchange : Pd adatom (fcc hollow)  ->  Pd swapped into the top layer + the displaced Au as adatom
             (same atoms both times: 45 Au + 1 Pd)
  depth    : Pd replacing one Au in layer 1 (surface), 2 (subsurface) or 3  (44 Au + 1 Pd)
usage:  python pd_alloy.py lj|emt|eam|mace|mpa|chgnet      ->  alloy_<pot>.csv
"""
import sys, csv, warnings
import numpy as np
from ase.build import bulk, fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.eos import EquationOfState
from ase.optimize import BFGS
from potentials import get_calc
warnings.filterwarnings("ignore")

pot = sys.argv[1] if len(sys.argv) > 1 else "emt"
calc = get_calc(pot)
vols, ens = [], []
for a in np.linspace(3.90, 4.30, 9):
    au = bulk("Au", "fcc", a=a, cubic=True); au.calc = calc
    vols.append(au.get_volume()); ens.append(au.get_potential_energy())
a0 = EquationOfState(vols, ens).fit()[0] ** (1 / 3)


def relax(atoms):
    atoms.set_constraint(FixAtoms(mask=[t > 3 for t in atoms.get_tags()]))  # freeze bottom 2 of 5
    atoms.calc = calc
    BFGS(atoms, logfile=None).run(fmax=0.02, steps=800)
    return atoms.get_potential_energy()


def central(slab, layer):
    """index of the atom in `layer` (tag) closest to the cell centre in x, y"""
    c = slab.cell[:2, :2].sum(axis=0) / 2
    idx = [i for i in range(len(slab)) if slab[i].tag == layer]
    return min(idx, key=lambda i: np.linalg.norm(slab.positions[i, :2] - c))


E = {}
# 1. adatom vs exchange
base = fcc111("Au", size=(3, 3, 5), a=a0, vacuum=12.0)
s = base.copy(); add_adsorbate(s, "Pd", 2.0, "fcc"); E["adatom"] = relax(s)
s = base.copy(); add_adsorbate(s, "Au", 2.0, "fcc")
s[central(s, 1)].symbol = "Pd"; E["exchange"] = relax(s)
# 2. substitution depth profile (45 atoms: 44 Au + 1 Pd)
for L in (1, 2, 3):
    s = base.copy(); s[central(s, L)].symbol = "Pd"; E[f"L{L}"] = relax(s)

ex = E["exchange"] - E["adatom"]
seg2, seg3 = E["L2"] - E["L1"], E["L3"] - E["L1"]
print(f"[{pot}] a0 = {a0:.3f} A")
print(f"exchange   (Pd into top layer, Au pops out)  dE = {ex:+.3f} eV  -> {'Pd sinks in' if ex < 0 else 'Pd stays on top'}")
print(f"subsurface (Pd layer 2 vs layer 1)           dE = {seg2:+.3f} eV  -> {'prefers subsurface' if seg2 < 0 else 'prefers surface'}")
print(f"layer 3    (Pd layer 3 vs layer 1)           dE = {seg3:+.3f} eV")
with open(f"alloy_{pot}.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["potential", "dE_exchange_eV", "dE_L2_minus_L1_eV", "dE_L3_minus_L1_eV"])
    w.writerow([pot, round(ex, 4), round(seg2, 4), round(seg3, 4)])
