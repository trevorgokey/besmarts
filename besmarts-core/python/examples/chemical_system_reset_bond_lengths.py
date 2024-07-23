"""
examples/chemical_system_reset_bond_lengths
"""


from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import smirnoff_models
from besmarts.perception import perception_rdkit
import pprint

pcp = perception_rdkit.perception_model_rdkit()
csys = smirnoff_models.smirnoff_load("./openff-2.1.0.offxml", pcp)

pos, extras = pcp.gcd.sdf_decode("./mol.sdf")
psys = mm.chemical_system_to_physical_system(csys, [pos])

bonds = mm.chemical_system_reset_bond_lengths(csys, [psys])
mm.chemical_system_print(csys, show_parameters=set((k[2] for k in bonds)))

