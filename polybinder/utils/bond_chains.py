import gsd.hoomd
from cmeutils.gsd_utils import get_molecule_cluster
import numpy as np


def bond_self_chains(snapshot, head_index, tail_index):
    cluster, clprops = get_molecule_cluster(snap=snapshot)
    num_new_bonds = len(cluster.cluster_keys)
    snapshot.bonds.N += num_new_bonds
    # Add to bonds group and ids
    new_groups = np.array(
        [np.array([key[0], key[-1]]) for key in cluster.cluster_keys]
    )
    new_ids = np.zeros(num_new_bonds)

    new_bond_groups = np.concatenate((snapshot.bonds.group, new_groups))
    new_bond_ids = np.concatenate((snapshot.bonds.typeid, new_ids))

    snapshot.bonds.group = new_bond_groups
    snapshot.bonds.typeid = new_bond_ids
    
    with gsd.hoomd.open("self-bonded.gsd", "wb") as new_traj:
        new_traj.append(snapshot)
    return snapshot

