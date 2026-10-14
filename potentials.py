"""The ladder of interatomic potentials for Pd on gold, all as ASE calculators.
    lj      12-6 Lennard-Jones, Heinz et al. 2008 metal parameters, Lorentz-Berthelot mixing
    emt     ASE's built-in Effective Medium Theory (Au and Pd are both native)
    eam     Embedded-Atom Method, Zhou-Johnson-Wadley 2004 Au-Pd alloy file (AuPd.eam.alloy)
    mace    MACE-MP-0 (medium)           universal machine-learned potential
    mpa     MACE-MPA-0 (medium)          same architecture, trained on 3x more data
    chgnet  CHGNet                       universal machine-learned potential
"""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import neighbor_list

KCAL = 0.0433641  # kcal/mol -> eV
# H. Heinz, R. A. Vaia, B. L. Farmer, R. R. Naik, J. Phys. Chem. C 112, 17281 (2008)
# Heinz writes E = eps[(r0/r)^12 - 2(r0/r)^6]; the usual 4*eps form needs sigma = r0 / 2^(1/6)
HEINZ = {"Au": (2.951 / 2 ** (1 / 6), 5.29 * KCAL), "Pd": (2.819 / 2 ** (1 / 6), 6.15 * KCAL)}  # sigma (A), eps (eV)


class MixedLJ(Calculator):
    """12-6 Lennard-Jones with per-element sigma/epsilon (ASE's LennardJones has only one pair)."""
    implemented_properties = ["energy", "forces"]

    def __init__(self, params=HEINZ, rc=10.0, **kw):
        super().__init__(**kw)
        self.params, self.rc = params, rc

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        i, j, d, D = neighbor_list("ijdD", self.atoms, self.rc)
        sym = self.atoms.get_chemical_symbols()
        s = np.array([(self.params[sym[a]][0] + self.params[sym[b]][0]) / 2 for a, b in zip(i, j)])
        e = np.array([np.sqrt(self.params[sym[a]][1] * self.params[sym[b]][1]) for a, b in zip(i, j)])
        if len(i) == 0:
            self.results = {"energy": 0.0, "forces": np.zeros((len(self.atoms), 3))}
            return
        sr6 = (s / d) ** 6
        sc6 = (s / self.rc) ** 6
        pair = 4 * e * (sr6**2 - sr6) - 4 * e * (sc6**2 - sc6)       # shifted to 0 at rc
        dEdr = 4 * e * (-12 * sr6**2 + 6 * sr6) / d
        f = (dEdr / d)[:, None] * D                                     # force on atom i
        forces = np.zeros((len(self.atoms), 3))
        np.add.at(forces, i, f)
        self.results = {"energy": 0.5 * pair.sum(), "forces": forces}


def get_calc(name):
    if name == "lj":
        return MixedLJ()
    if name == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    if name == "eam":
        from ase.calculators.eam import EAM
        return EAM(potential="AuPd.eam.alloy")
    if name in ("mace", "mpa"):
        from mace.calculators import mace_mp
        model = "medium" if name == "mace" else "medium-mpa-0"
        return mace_mp(model=model, default_dtype="float64", device="cpu")
    if name == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        return CHGNetCalculator(use_device="cpu")
    raise SystemExit("potential must be one of: lj emt eam mace mpa chgnet")
