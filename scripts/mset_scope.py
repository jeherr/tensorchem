import glob
from tensorchem.dataset.molecule import MoleculeSet as MSet

def get_mset_info(mset_mols):
    msets = []
    mset_natoms = []
    mset_elements = []
    for mol in mset_mols:
        mset = MSet()
        mset.filename = mol
        mset.load()

        repeat = 0
        #Removing Non-Unique Molecules
        for prev_mset in msets:
            if mset.is_isomer(prev_mset):
                repeat = 1
        
        #Removing Metals
        try:
            for element_set in mset.elements:
                for element in element_set:
                    if element=='Li' or element=='Be' or element=='As' or element=='K' or element=='Na' or element=='Ca' or element=='Mg':
                        repeat = 1
        except:
            pass

        if repeat == 0:
            msets.append(mset)
            mset_natoms.append(mset.n_atoms)
            mset_elements.append(mset.elements)

    elements = {'H': 0, 'C': 0, 'B': 0, 'N': 0, 'O': 0, 'F': 0, 'Si': 0, 'P': 0, 'S': 0, 'Cl': 0, 'Se': 0, 'Br': 0, 'I': 0}
    for mol_elements in mset_elements:
        for element_set in mol_elements:
            for element in element_set:
                elements[element] += 1

    for element, value in elements.items():
        if value == 0:
            elements.pop(element)

    return (msets, mset_natoms, elements)

def get_meta_mset_info():
    meta_msets = []
    meta_mset_natoms = []
    meta_mset_elements = []
    for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta/*.mset")+glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/aimd/*.mset"):
        mset = MSet()
        mset.filename = mol
        mset.load()

        repeat = 0
        #Removing Non-Unique Molecules
        for prev_mset in meta_msets:
            if mset.is_isomer(prev_mset):
                repeat = 1
        
        #Removing Metals
        for element_set in mset.elements:
            for element in element_set:
                if element=='Li' or element=='Be' or element=='As' or element=='K' or element=='Na' or element=='Ca' or element=='Mg':
                    repeat = 1

        if repeat == 0:
            meta_msets.append(mset)
            meta_mset_natoms.append(mset.n_atoms)
            meta_mset_elements.append(mset.elements)

    meta_elements = {'H': 0, 'C': 0, 'B': 0, 'N': 0, 'O': 0, 'F': 0, 'Si': 0, 'P': 0, 'S': 0, 'Cl': 0, 'Se': 0, 'Br': 0, 'I': 0}
    for mol_elements in meta_mset_elements:
        for element_set in mol_elements:
            for element in element_set:
                meta_elements[element] += 1

    for element, value in meta_elements.items():
        if value == 0:
            meta_elements.pop(element)

    return (meta_msets, meta_mset_natoms, meta_elements)

def get_opt_mset_info():
    opt_msets = []
    opt_mset_natoms = []
    opt_mset_elements = []
    for i in range(1,10):
        for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt/"+str(i)+"/*.mset"):
            mset = MSet()
            mset.filename = mol
            mset.load()

            repeat = 0
            #Removing Non-Unique Molecules
            for prev_mset in opt_msets:
                if mset.is_isomer(prev_mset):
                    repeat = 1
            try:
                for element_set in mset.elements:
                    for element in element_set:
                        if element=='Li' or element=='Be' or element=='As' or element=='K' or element=='Na' or element=='Ca' or element=='Mg':
                            repeat = 1
            except:
                pass
        
            if repeat == 0:
                opt_msets.append(mset)
                opt_mset_natoms.append(mset.n_atoms)
                opt_mset_elements.append(mset.elements)

    for x in ['b', 'br', 'i', 'p', 'p_new', 'se', 'si']:
        for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_cs40/"+str(x)+"/*.mset"):
            mset = MSet()
            mset.filename = mol
            mset.load()
        
            try:
                for element_set in mset.elements:
                    for element in element_set:
                        if element=='Li' or element=='Be' or element=='As' or element=='K' or element=='Na' or element=='Ca' or element=='Mg':
                            repeat = 1
            except:
                pass

            repeat = 0
            #Removing Non-Unique Molecules
            for prev_mset in opt_msets:
                if mset.is_isomer(prev_mset):
                    repeat = 1
        
            if repeat == 0:
                opt_msets.append(mset)
                opt_mset_natoms.append(mset.n_atoms)
                opt_mset_elements.append(mset.elements)

    opt_elements = {'H': 0, 'C': 0, 'B': 0, 'N': 0, 'O': 0, 'F': 0, 'Si': 0, 'P': 0, 'S': 0, 'Cl': 0, 'Se': 0, 'Br': 0, 'I': 0}
    for mol_elements in opt_mset_elements:
        for element_set in mol_elements:
            for element in element_set:
                opt_elements[element] += 1

    for element, value in opt_elements.items():
        if value == 0:
            opt_elements.pop(element)

    return (opt_msets, opt_mset_natoms, opt_elements)

def get_ani_mset_info():
    ani_msets = []
    ani_mset_natoms = []
    ani_mset_elements = []

    for mol in glob.glob("/mnt/sdb1/adriscoll/ani1x-data/ani1x-msets/*.mset"):
        mset = MSet()
        mset.filename = mol
        mset.load()

        repeat = 0
        #Removing Non-Unique Molecules
        for prev_mset in ani_msets:
            if mset.is_isomer(prev_mset):
                 repeat = 1
        
        if repeat == 0:
            ani_msets.append(mset)
            ani_mset_natoms.append(mset.n_atoms)
            ani_mset_elements.append(mset.elements)

    ani_elements = {'H': 0, 'C': 0, 'N': 0, 'O': 0}
    for mol_elements in ani_mset_elements:
        for element_set in mol_elements:
            for element in element_set:
                ani_elements[element] += 1

    for element, value in ani_elements.items():
        if value == 0:
            ani_elements.pop(element)

    return (ani_msets, ani_mset_natoms, ani_elements)
