"""
This is shamelessly taken from https://github.com/jensengroup/xyz2mol with credit to the authors of the work at
DOI: 10.1002/bkcs.10334 as well for the algorithm. This is a sometimes useful tool to try and generate and RDKit
Mol object from a set of XYZ coordinates of a molecule.
"""

import copy
import itertools

try:
    from rdkit.Chem import rdEHTTools  # requires RDKit 2019.9.1 or later
except ImportError:
    rdEHTTools = None

from collections import defaultdict

import numpy as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem


def get_UA(max_valence_list, valence_list):
    """
    """
    UA = []
    DU = []
    for i, (max_valence, valence) in enumerate(zip(max_valence_list, valence_list)):
        if not max_valence - valence > 0:
            continue
        UA.append(i)
        DU.append(max_valence - valence)
    return UA, DU


def get_BO(AC, DU, valences, UA_pairs, use_graph=True):
    """
    """
    BO = AC.copy()
    DU_save = []

    while DU_save != DU:
        for i, j in UA_pairs:
            BO[i, j] += 1
            BO[j, i] += 1

        BO_valence = list(BO.sum(axis=1))
        DU_save = copy.copy(DU)
        UA, DU = get_UA(valences, BO_valence)
        UA_pairs = get_UA_pairs(UA, AC, use_graph=use_graph)[0]

    return BO


def valences_not_too_large(BO, valences):
    """
    """
    number_of_bonds_list = BO.sum(axis=1)
    for valence, number_of_bonds in zip(valences, number_of_bonds_list):
        if number_of_bonds > valence:
            return False

    return True


def charge_is_ok(BO, charge, atomic_valence_electrons, atoms, allow_charged_fragments=True):
    # total charge
    Q = 0

    # charge fragment list
    q_list = []

    if allow_charged_fragments:

        BO_valences = list(BO.sum(axis=1))
        for i, atom in enumerate(atoms):
            q = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])
            Q += q
            if atom == 6:
                num_single_bonds_to_carbon = list(BO[i, :]).count(1)
                if num_single_bonds_to_carbon == 2 and BO_valences[i] == 2:
                    Q += 1
                    q = 2
                if num_single_bonds_to_carbon == 3 and Q + 1 < charge:
                    Q += 2
                    q = 1

            if q != 0:
                q_list.append(q)

    return charge == Q


def BO_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
             allow_charged_fragments=True):
    """
    Sanity of bond-orders
    args:
        BO -
        AC -
        charge -
        DU -
    optional
        allow_charges_fragments -
    returns:
        boolean - true of molecule is OK, false if not
    """

    if not valences_not_too_large(BO, valences):
        return False

    check_sum = (BO - AC).sum() == sum(DU)
    check_charge = charge_is_ok(BO, charge, atomic_valence_electrons, atoms, allow_charged_fragments)

    if check_charge and check_sum:
        return True

    return False


def get_atomic_charge(atom, atomic_valence_electrons, BO_valence):
    """
    """

    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence

    return charge


def clean_charges(mol):
    """
    This hack should not be needed anymore, but is kept just in case
    """

    Chem.SanitizeMol(mol)
    # rxn_smarts = ['[N+:1]=[*:2]-[C-:3]>>[N+0:1]-[*:2]=[C-0:3]',
    #              '[N+:1]=[*:2]-[O-:3]>>[N+0:1]-[*:2]=[O-0:3]',
    #              '[N+:1]=[*:2]-[*:3]=[*:4]-[O-:5]>>[N+0:1]-[*:2]=[*:3]-[*:4]=[O-0:5]',
    #              '[#8:1]=[#6:2]([!-:6])[*:3]=[*:4][#6-:5]>>[*-:1][*:2]([*:6])=[*:3][*:4]=[*+0:5]',
    #              '[O:1]=[c:2][c-:3]>>[*-:1][*:2][*+0:3]',
    #              '[O:1]=[C:2][C-:3]>>[*-:1][*:2]=[*+0:3]']

    rxn_smarts = [''.join(('[#6,#7:1]1=[#6,#7:2][#6,#7:3]=[#6,#7:4][CX3-,NX3-:5][#6,#7:6]1=[#6,#7:7]>>',
                           '[#6,#7:1]1=[#6,#7:2][#6,#7:3]=[#6,#7:4][-0,-0:5]=[#6,#7:6]1[#6-,#7-:7]')),
                  ''.join(('[#6,#7:1]1=[#6,#7:2][#6,#7:3](=[#6,#7:4])[#6,#7:5]=[#6,#7:6][CX3-,NX3-:7]1>>',
                           '[#6,#7:1]1=[#6,#7:2][#6,#7:3]([#6-,#7-:4])=[#6,#7:5][#6,#7:6]=[-0,-0:7]1'))]

    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    for i, fragment in enumerate(fragments):
        for smarts in rxn_smarts:
            patt = Chem.MolFromSmarts(smarts.split(">>")[0])
            while fragment.HasSubstructMatch(patt):
                rxn = AllChem.ReactionFromSmarts(smarts)
                ps = rxn.RunReactants((fragment,))
                fragment = ps[0][0]
                Chem.SanitizeMol(fragment)
        if i == 0:
            mol = fragment
        else:
            mol = Chem.CombineMols(mol, fragment)

    return mol


def BO2mol(mol, BO_matrix, atoms, atomic_valence_electrons,
           mol_charge, allow_charged_fragments=True):
    """
    based on code written by Paolo Toscani
    From bond order, atoms, valence structure and total charge, generate an
    rdkit molecule.
    args:
        mol - rdkit molecule
        BO_matrix - bond order matrix of molecule
        atoms - list of integer atomic symbols
        atomic_valence_electrons -
        mol_charge - total charge of molecule
    optional:
        allow_charged_fragments - bool - allow charged fragments
    returns
        mol - updated rdkit molecule with bond connectivity
    """

    l = len(BO_matrix)
    l2 = len(atoms)
    BO_valences = list(BO_matrix.sum(axis=1))

    if l != l2:
        raise RuntimeError('sizes of adjMat ({0:d}) and Atoms {1:d} differ'.format(l, l2))

    rw_mol = Chem.RWMol(mol)

    bond_type_dict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }

    for i in range(l):
        for j in range(i + 1, l):
            bo = int(round(BO_matrix[i, j]))
            if bo == 0:
                continue
            bt = bond_type_dict.get(bo, Chem.BondType.SINGLE)
            rw_mol.AddBond(i, j, bt)

    mol = rw_mol.GetMol()

    if allow_charged_fragments:
        mol = set_atomic_charges(
            mol,
            atoms,
            atomic_valence_electrons,
            BO_valences,
            BO_matrix,
            mol_charge)
    else:
        mol = set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences)

    return mol


def set_atomic_charges(mol, atoms, atomic_valence_electrons,
                       BO_valences, BO_matrix, mol_charge):
    """
    """
    q = 0
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])
        q += charge
        if atom == 6:
            num_single_bonds_to_carbon = list(BO_matrix[i, :]).count(1)
            if num_single_bonds_to_carbon == 2 and BO_valences[i] == 2:
                q += 1
                charge = 0
            if num_single_bonds_to_carbon == 3 and q + 1 < mol_charge:
                q += 2
                charge = 1

        if abs(charge) > 0:
            a.SetFormalCharge(int(charge))

    mol = clean_charges(mol)

    return mol


def set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences):
    """
    The number of radical electrons = absolute atomic charge
    """
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(
            atom,
            atomic_valence_electrons[atom],
            BO_valences[i])

        if abs(charge) > 0:
            a.SetNumRadicalElectrons(abs(int(charge)))

    return mol


def get_bonds(UA, AC):
    """
    """
    bonds = []

    for k, i in enumerate(UA):
        for j in UA[k + 1:]:
            if AC[i, j] == 1:
                bonds.append(tuple(sorted([i, j])))

    return bonds


def get_UA_pairs(UA, AC, use_graph=True):
    """
    """

    bonds = get_bonds(UA, AC)

    if len(bonds) == 0:
        return [()]

    if use_graph:
        G = nx.Graph()
        G.add_edges_from(bonds)
        UA_pairs = [list(nx.max_weight_matching(G))]
        return UA_pairs

    max_atoms_in_combo = 0
    UA_pairs = [()]
    for combo in list(itertools.combinations(bonds, int(len(UA) / 2))):
        flat_list = [item for sublist in combo for item in sublist]
        atoms_in_combo = len(set(flat_list))
        if atoms_in_combo > max_atoms_in_combo:
            max_atoms_in_combo = atoms_in_combo
            UA_pairs = [combo]

        elif atoms_in_combo == max_atoms_in_combo:
            UA_pairs.append(combo)

    return UA_pairs


def AC2BO(AC, atoms, charge, allow_charged_fragments=True, use_graph=True):
    """
    implemenation of algorithm shown in Figure 2
    UA: unsaturated atoms
    DU: degree of unsaturation (u matrix in Figure)
    best_BO: Bcurr in Figure
    """
    atomic_valence = defaultdict(list)
    atomic_valence[1] = [1]
    atomic_valence[5] = [3, 4]
    atomic_valence[6] = [4]
    atomic_valence[7] = [3, 4]
    atomic_valence[8] = [2, 1, 3]
    atomic_valence[9] = [1]
    atomic_valence[14] = [4]
    atomic_valence[15] = [5, 3]  # [5,4,3]
    atomic_valence[16] = [6, 3, 2]  # [6,4,2]
    atomic_valence[17] = [1]
    atomic_valence[32] = [4]
    atomic_valence[35] = [1]
    atomic_valence[53] = [1]

    atomic_valence_electrons = {}
    atomic_valence_electrons[1] = 1
    atomic_valence_electrons[2] = 2
    atomic_valence_electrons[3] = 1
    atomic_valence_electrons[4] = 2
    atomic_valence_electrons[5] = 3
    atomic_valence_electrons[6] = 4
    atomic_valence_electrons[7] = 5
    atomic_valence_electrons[8] = 6
    atomic_valence_electrons[9] = 7
    atomic_valence_electrons[11] = 1
    atomic_valence_electrons[12] = 2
    atomic_valence_electrons[13] = 2
    atomic_valence_electrons[14] = 4
    atomic_valence_electrons[15] = 5
    atomic_valence_electrons[16] = 6
    atomic_valence_electrons[17] = 7
    atomic_valence_electrons[19] = 1
    atomic_valence_electrons[20] = 2
    atomic_valence_electrons[31] = 3
    atomic_valence_electrons[32] = 4
    atomic_valence_electrons[33] = 5
    atomic_valence_electrons[34] = 6
    atomic_valence_electrons[35] = 7
    atomic_valence_electrons[37] = 1
    atomic_valence_electrons[38] = 2
    atomic_valence_electrons[49] = 3
    atomic_valence_electrons[50] = 4
    atomic_valence_electrons[51] = 5
    atomic_valence_electrons[52] = 6
    atomic_valence_electrons[53] = 7

    # make a list of valences, e.g. for CO: [[4],[2,1]]
    valences_list_of_lists = []
    AC_valence = list(AC.sum(axis=1))

    for atomicNum, valence in zip(atoms, AC_valence):
        # valence can't be smaller than number of neighbourgs
        possible_valence = [x for x in atomic_valence[atomicNum] if x >= valence]
        valences_list_of_lists.append(possible_valence)

    # convert [[4],[2,1]] to [[4,2],[4,1]]
    valences_list = itertools.product(*valences_list_of_lists)

    best_BO = AC.copy()
    charge_ok = False

    for valences in valences_list:

        UA, DU_from_AC = get_UA(valences, AC_valence)

        check_len = (len(UA) == 0)
        if check_len:
            check_bo = BO_is_OK(AC, AC, charge, DU_from_AC,
                                atomic_valence_electrons, atoms, valences,
                                allow_charged_fragments=allow_charged_fragments)
        else:
            check_bo = None

        if check_len and check_bo:
            return AC, atomic_valence_electrons

        UA_pairs_list = get_UA_pairs(UA, AC, use_graph=use_graph)
        for UA_pairs in UA_pairs_list:
            BO = get_BO(AC, DU_from_AC, valences, UA_pairs, use_graph=use_graph)
            status = BO_is_OK(BO, AC, charge, DU_from_AC,
                              atomic_valence_electrons, atoms, valences,
                              allow_charged_fragments=allow_charged_fragments)
            charge_ok = charge_is_ok(BO, charge, atomic_valence_electrons, atoms,
                                     allow_charged_fragments=allow_charged_fragments)

            if status:
                return BO, atomic_valence_electrons
            elif BO.sum() >= best_BO.sum() and valences_not_too_large(BO, valences) and charge_ok:
                best_BO = BO.copy()

    if not charge_ok:
        print("Warning: SMILES charge doesn't match input charge")
    return best_BO, atomic_valence_electrons


def AC2mol(mol, AC, atoms, charge, allow_charged_fragments=True, use_graph=True):
    """
    """

    # convert AC matrix to bond order (BO) matrix
    BO, atomic_valence_electrons = AC2BO(
        AC,
        atoms,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph)

    # add BO connectivity and charge info to mol object
    mol = BO2mol(
        mol,
        BO,
        atoms,
        atomic_valence_electrons,
        charge,
        allow_charged_fragments=allow_charged_fragments)

    return mol


def get_proto_mol(atoms):
    """
    """
    mol = Chem.MolFromSmarts("[#" + str(atoms[0]) + "]")
    rwMol = Chem.RWMol(mol)
    for i in range(1, len(atoms)):
        a = Chem.Atom(atoms[i])
        rwMol.AddAtom(a)

    mol = rwMol.GetMol()

    return mol


def xyz_to_connectivity(atoms, xyz, charge, use_huckel=False):
    """
    atoms and coordinates to atom connectivity (AC)
    args:
        atoms - int atom types
        xyz - coordinates
        charge - molecule charge
    optional:
        use_huckel - Use Huckel method for atom connecitivty
    returns
        ac - atom connectivity matrix
        mol - rdkit molecule
    """

    if use_huckel:
        return xyz2AC_huckel(atoms, xyz, charge)
    else:
        return xyz2AC_vdW(atoms, xyz)


def xyz2AC_vdW(atoms, xyz):
    # Get mol template
    mol = get_proto_mol(atoms)

    # Set coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)

    AC = get_AC(mol)

    return AC, mol


def get_AC(mol, covalent_factor=1.3):
    """
    Generate adjacent matrix from atoms and coordinates.
    AC is a (num_atoms, num_atoms) matrix with 1 being covalent bond and 0 is not
    covalent_factor - 1.3 is an arbitrary factor
    args:
        mol - rdkit molobj with 3D conformer
    optional
        covalent_factor - increase covalent bond length threshold with facto
    returns:
        AC - adjacent matrix
    """

    # Calculate distance matrix
    dMat = Chem.Get3DDistanceMatrix(mol)

    pt = Chem.GetPeriodicTable()
    num_atoms = mol.GetNumAtoms()
    AC = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum()) * covalent_factor
        for j in range(i + 1, num_atoms):
            a_j = mol.GetAtomWithIdx(j)
            Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum()) * covalent_factor
            if dMat[i, j] <= Rcov_i + Rcov_j:
                AC[i, j] = 1
                AC[j, i] = 1

    return AC


def xyz2AC_huckel(atomicNumList, xyz, charge):
    """
    args
        atomicNumList - atom type list
        xyz - coordinates
        charge - molecule charge
    returns
        ac - atom connectivity
        mol - rdkit molecule
    """
    mol = get_proto_mol(atomicNumList)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)

    num_atoms = len(atomicNumList)
    AC = np.zeros((num_atoms, num_atoms)).astype(int)

    mol_huckel = Chem.Mol(mol)
    mol_huckel.GetAtomWithIdx(0).SetFormalCharge(charge)  # mol charge arbitrarily added to 1st atom

    passed, result = rdEHTTools.RunMol(mol_huckel)
    opop = result.GetReducedOverlapPopulationMatrix()
    tri = np.zeros((num_atoms, num_atoms))
    tri[np.tril(np.ones((num_atoms, num_atoms), dtype=bool))] = opop  # lower triangular to square matrix
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            pair_pop = abs(tri[j, i])
            if pair_pop >= 0.15:  # arbitry cutoff for bond. May need adjustment
                AC[i, j] = 1
                AC[j, i] = 1

    return AC, mol


def chiral_stereo_check(mol):
    """
    Find and embed chiral information into the model based on the coordinates
    args:
        mol - rdkit molecule, with embeded conformer
    """
    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol, -1)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(mol, -1)

    return


def geom_to_mol(geom, charge, allow_charged_fragments=True, use_graph=True, use_huckel=False, embed_chiral=True):
    """
    Generate a rdkit molobj from atoms, coordinates and a total_charge.
    args:
        atoms - list of atom types (int)
        coordinates - 3xN Cartesian coordinates
        charge - total charge of the system (default: 0)
    optional:
        allow_charged_fragments - alternatively radicals are made
        use_graph - use graph (networkx)
        use_huckel - Use Huckel method for atom connectivity prediction
        embed_chiral - embed chiral information to the molecule
    returns:
        mol - rdkit molobj
    """

    # Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge,
    # and mol object with no connectivity information
    atoms = [atom.at_num for atom in geom.atoms]
    AC, mol = xyz_to_connectivity(atoms, geom.coords, charge, use_huckel=use_huckel)

    # Convert AC to bond order matrix and add connectivity and charge info to
    # mol object
    new_mol = AC2mol(mol, AC, atoms, charge,
                     allow_charged_fragments=allow_charged_fragments,
                     use_graph=use_graph)

    # Check for stereocenters and chiral centers
    if embed_chiral:
        chiral_stereo_check(new_mol)

    return new_mol


def geom_to_smiles(geom, charge=0):
    m = geom_to_mol(geom, charge)
    smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    m = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    return smiles
