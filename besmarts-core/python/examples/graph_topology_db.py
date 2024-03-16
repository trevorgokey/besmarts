"""
"""

from besmarts.core import assignments
from besmarts.core import graphs
from besmarts.data import smiles_graphs

pos = smiles_graphs.ethane_pos()
db = assignments.graph_topology_db()
gid = assignments.graph_topology_db_add_graph(db, graphs.graph({}, {}))
assignments.graph_topology_db_add_positions(db, gid, pos.selections)

sid = (1,)
db.gfuncs[POSITIONS].selection[gid][sid]
db.gfuncs[PHYSICAL_MODEL_BONDS_LABELS].selection[gid][sid][tid][pid] 
db.gfuncs[PHYSICAL_MODEL_BONDS_VALUES].selection[gid][sid][tid][pid]

db.gfuncs[PHYSICAL_MODEL_BONDS_VALUES].selection[gid][sid] = {tid: pid}


