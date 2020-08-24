import glob
import os
from ase.data import atomic_numbers
from tensorchem.dataset.molecule import MoleculeSet, Geometry, Atom


def parse_aimd_input(f):
    n_atoms = 0
    mol_flag = True
    job_flag = True
    while mol_flag or job_flag:
        line = next(f)
        if "$molecule" in line:
            line = next(f)
            spin, multiplicity = int(line.split()[0]), int(line.split()[1])
            while True:
                line = next(f)
                if "$end" in line:
                    mol_flag = False
                    break
                if line:
                    at_sym = line.split()[0]
                    n_atoms += 1
        elif "$rem" in line:
            while True:
                line = next(f)
                if "$end" in line:
                    job_flag = False
                    break
                elif "method" in line:
                    method = line.split()[-1]
                elif "basis" in line:
                    basis = line.split()[-1]
                elif "time_step" in line:
                    time_step = int(line.split()[-1])
                elif "aimd_steps" in line:
                    aimd_steps = int(line.split()[-1])
                elif "aimd_temp" in line:
                    aimd_temp = int(line.split()[1])
    return n_atoms, at_sym, method, basis, time_step, aimd_steps, aimd_temp


def parse_sp_input(f):
    n_atoms = 0
    mol_flag = True
    job_flag = True
    while mol_flag or job_flag:
        line = next(f)
        if "$molecule" in line:
            line = next(f)
            spin, multiplicity = int(line.split()[0]), int(line.split()[1])
            while True:
                line = next(f)
                if "$end" in line:
                    mol_flag = False
                    break
                if line:
                    at_sym = line.split()[0]
                    n_atoms += 1
        elif "$rem" in line:
            while True:
                line = next(f)
                if "$end" in line:
                    job_flag = False
                    break
                elif "method" in line:
                    method = line.split()[-1]
                elif "basis" in line:
                    basis = line.split()[-1]
    return n_atoms, at_sym, method, basis


def parse_energy(f):
    while True:
        line = next(f)
        if "Convergence criterion met" in line:
            break
    return float(line.split()[1])


def parse_forces(f, n_atoms):
    force = []
    while len(force) < n_atoms:
        line = next(f)
        new_atoms = len(line.split())
        new_forces = [[] for _ in range(new_atoms)]
        for _ in range(3):
            line = next(f).split()
            for atm, frc in enumerate(line[1:]):
                try:
                    new_forces[atm].append(float(frc) / -0.529177208590000)
                except ValueError as ex:
                    print(ex)
                    print("Probably force gradient was to large to be parsed correctly, skipping this geometry")
                    return None
        force += new_forces
    return force


def parse_atoms_coords(f, n_atoms):
    next(f)
    next(f)
    atom_nums, coord = [], []
    for _ in range(n_atoms):
        split_line = next(f).split()
        at_sym = split_line[1]
        at_num = atomic_numbers[at_sym]
        atom_nums.append(at_num)
        try:
            coord.append((float(split_line[2]), float(split_line[3]), float(split_line[4])))
        except ValueError as ex:
            print(ex)
            print("Weird starting coordinates probably, skipping this geometry")
            return None, None
    return atom_nums, coord


def parse_dipole(f):
    line = next(f).split()
    dipole = [float(comp) for comp in line[1::2]]
    return dipole


def parse_quadrupole(f):
    quadrupole = []
    for _ in range(2):
        line = next(f).split()
        quadrupole += [float(comp) for comp in line[1::2]]
    return quadrupole


def parse_charges(f, n_atoms):
    charge = []
    for _ in range(3):
        next(f)
    for _ in range(n_atoms):
        line = next(f).split()
        charge.append(float(line[-1]))
    return charge


def build_new_geom(atomic_nums, coords, energy, forces, dipole, quadrupole, charges, method, basis):
    mol_labels = {".".join(("potential", method, basis)): energy[-1],
                  ".".join(("dipole", method, basis)): dipole[-1],
                  ".".join(("quadrupole", method, basis)): quadrupole[-1]}
    atom_jsons = [{'atomic_num': at_num} for at_num in atomic_nums[-1]]
    for i, atom in enumerate(atom_jsons):
        atom['xyz'] = coords[-1][i]
        atom['labels'] = {
            ".".join(("charge", "mulliken", method, basis)): charges[-1][i],
            ".".join(("forces", method, basis)): forces[-1][i]
        }
    atoms = [Atom.from_json(at_json) for at_json in atom_jsons]
    new_geom = Geometry(atoms, mol_labels)
    return new_geom


def read_aimd_data(filename):
    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []
    print("Reading ", filename)
    with open(filename, "r") as f:
        try:
            while True:
                line = next(f)
                if "User input:" in line:
                    n_atoms, at_sym, method, basis, time_steps, aimd_steps, aimd_temp = parse_aimd_input(f)
                # These next sections **SHOULD** only grab the first instance of these properties (i.e. before the AIMD
                # section starts)
                elif "Standard Nuclear Orientation" in line:
                    atom_nums, coord = parse_atoms_coords(f, n_atoms)
                    atomic_nums.append(atom_nums)
                    coords.append(coord)
                # energy
                elif "Cycle" in line and "Energy" in line:
                    energies.append(parse_energy(f))
                # forces
                elif "Gradient of SCF Energy" in line:
                    forces.append(parse_forces(f, n_atoms))
                # dipoles
                elif "Dipole Moment (Debye)" in line:
                    dipoles.append(parse_dipole(f))
                # quadrupoles
                elif "Quadrupole Moments (Debye-Ang)" in line:
                    quadrupoles.append(parse_quadrupole(f))
                # charges
                elif "Ground-State Mulliken Net Atomic Charges" in line:
                    charges.append(parse_charges(f, n_atoms))
                elif "AB INITIO MOLECULAR DYNAMICS" in line:
                    atoms = [Atom(at_num) for at_num in atomic_nums[0]]
                    mset = MoleculeSet(atoms)
                    aimd_trajectory = [
                        build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles, charges,
                                       method, basis)]
                    while True:
                        line = next(f)
                        if "TIME STEPS COMPLETED" in line:
                            if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                                    quadrupoles) == len(charges):
                                if len(energies) > len(aimd_trajectory):
                                    aimd_trajectory.append(
                                        build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles,
                                                       charges, method, basis))
                            else:
                                print("Error reading AIMD trajectory at last step", filename,
                                      "Returning MSet before this step")
                            mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                            return mset
                        elif "TIME STEP" in line:
                            time_step = int(line.split()[2].lstrip("#"))
                            time_au = float(line.split()[5])
                            time_fs = float(line.split()[8])
                            if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                                    quadrupoles) == len(charges):
                                aimd_trajectory.append(
                                    build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles, charges,
                                                   method, basis))
                            else:
                                print("Error reading AIMD trajectory at step", time_step, filename,
                                      "Returning MSet before "
                                      "this step")
                                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                                return mset
                        elif "Standard Nuclear Orientation" in line:
                            atom_nums, coord = parse_atoms_coords(f, n_atoms)
                            atomic_nums.append(atom_nums)
                            coords.append(coord)
                        # energy
                        elif "Cycle" in line and "Energy" in line:
                            energies.append(parse_energy(f))
                        # forces
                        elif "Gradient of SCF Energy" in line:
                            forces.append(parse_forces(f, n_atoms))
                        # dipoles
                        elif "Dipole Moment (Debye)" in line:
                            dipoles.append(parse_dipole(f))
                        # quadrupoles
                        elif "Quadrupole Moments (Debye-Ang)" in line:
                            quadrupoles.append(parse_quadrupole(f))
                        # charges
                        elif "Ground-State Mulliken Net Atomic Charges" in line:
                            charges.append(parse_charges(f, n_atoms))
        except StopIteration:
            print("Hit EOF on ", filename)
            try:
                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                return mset
            except UnboundLocalError:
                print("No MSet built for ", filename)
                return None


def read_sp_data(filename):
    n_atoms = None
    mset = None
    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []
    print("Reading ", filename)
    with open(filename, "r") as f:
        try:
            while True:
                line = next(f)
                if "User input:" in line:
                    n_atoms, at_sym, method, basis = parse_sp_input(f)
                elif "Standard Nuclear Orientation" in line:
                    atom_nums, coord = parse_atoms_coords(f, n_atoms)
                    atomic_nums.append(atom_nums)
                    coords.append(coord)
                # energy
                elif "Cycle" in line and "Energy" in line:
                    energies.append(parse_energy(f))
                # forces
                elif "Gradient of SCF Energy" in line:
                    force = parse_forces(f, n_atoms)
                    if force is not None:
                        forces.append(force)
                    else:
                        print(filename, " contains unparsed forces")
                # dipoles
                elif "Dipole Moment (Debye)" in line:
                    dipoles.append(parse_dipole(f))
                # quadrupoles
                elif "Quadrupole Moments (Debye-Ang)" in line:
                    quadrupoles.append(parse_quadrupole(f))
                # charges
                elif "Ground-State Mulliken Net Atomic Charges" in line:
                    charges.append(parse_charges(f, n_atoms))
                elif "Thank you very much for using Q-Chem." in line:
                    if mset is None:
                        print("Attempting to make initial MSet for ", filename)
                        try:
                            atoms = [Atom(at_num) for at_num in atomic_nums[0]]
                            mset = MoleculeSet(atoms)
                            meta_trajectory = []
                        except:
                            print("Error making MSet for ", filename)
                            return None
                    if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                            quadrupoles) == len(charges):
                        if len(energies) > len(meta_trajectory):
                            meta_trajectory.append(
                                build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles, charges,
                                               method, basis))
        except StopIteration:
            print("Hit EOF on ", filename)
            try:
                mset.trajectories[".".join((method, basis, "meta"))] = meta_trajectory
                return mset
            except UnboundLocalError:
                print("No MSet built for ", filename)
                return None


def read_multi_sp_data(filenames):
    n_atoms = None
    mset = None
    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []

    for filename in filenames:
        with open(filename, "r") as f:
            try:
                while True:
                    line = next(f)
                    if "User input:" in line:
                        n_atoms, at_sym, method, basis = parse_sp_input(f)
                    elif "Standard Nuclear Orientation" in line:
                        atom_nums, coord = parse_atoms_coords(f, n_atoms)
                        atomic_nums.append(atom_nums)
                        coords.append(coord)
                    # energy
                    elif "Cycle" in line and "Energy" in line:
                        energies.append(parse_energy(f))
                    # forces
                    elif "Gradient of SCF Energy" in line:
                        force = parse_forces(f, n_atoms)
                        if force is not None:
                            forces.append(force)
                        else:
                            print(filename, " contains unparsed forces")
                    # dipoles
                    elif "Dipole Moment (Debye)" in line:
                        dipoles.append(parse_dipole(f))
                    # quadrupoles
                    elif "Quadrupole Moments (Debye-Ang)" in line:
                        quadrupoles.append(parse_quadrupole(f))
                    # charges
                    elif "Ground-State Mulliken Net Atomic Charges" in line:
                        charges.append(parse_charges(f, n_atoms))
                    elif "Thank you very much for using Q-Chem." in line:
                        if mset is None:
                            print("Attempting to make initial MSet for ", filename)
                            try:
                                atoms = [Atom(at_num) for at_num in atomic_nums[0]]
                                mset = MoleculeSet(atoms)
                                meta_trajectory = []
                            except:
                                print("Error making MSet for ", filename)
                                return None
                        if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                                quadrupoles) == len(charges):
                            if len(energies) > len(meta_trajectory):
                                meta_trajectory.append(
                                    build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles, charges,
                                                   method, basis))
                        else:
                            min_len = min([len(atomic_nums), len(energies), len(forces), len(coords), len(dipoles), len(
                                quadrupoles), len(charges)])
                            atomic_nums, energies, forces, coords, dipoles, quadrupoles, charges = atomic_nums[
                                                                                                   :min_len], \
                                                                                                   energies[
                                                                                                   :min_len], forces[
                                                                                                              :min_len], coords[
                                                                                                                         :min_len], dipoles[
                                                                                                                                    :min_len], quadrupoles[
                                                                                                                                               :min_len], charges[
                                                                                                                                                          :min_len]
            except (StopIteration, UnicodeDecodeError):
                if mset is not None:
                    if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                            quadrupoles) == len(charges):
                        if len(energies) > len(meta_trajectory):
                            meta_trajectory.append(
                                build_new_geom(atomic_nums, coords, energies, forces, dipoles, quadrupoles, charges,
                                               method, basis))
                    else:
                        min_len = min([len(atomic_nums), len(energies), len(forces), len(coords), len(dipoles), len(
                            quadrupoles), len(charges)])
                        atomic_nums, energies, forces, coords, dipoles, quadrupoles, charges = atomic_nums[:min_len], \
                                                                                               energies[
                                                                                               :min_len], forces[
                                                                                                          :min_len], coords[
                                                                                                                     :min_len], dipoles[
                                                                                                                                :min_len], quadrupoles[
                                                                                                                                           :min_len], charges[
                                                                                                                                                      :min_len]
                continue
    if mset is None:
        print("No geometries collected for ", filenames[0])
        return None
    mset.trajectories[".".join((method, basis, "meta"))] = meta_trajectory
    return mset


def read_opt_data(filename):
    n_atoms = None
    mset = None
    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []
    print("Reading ", filename)
    with open(filename, "r") as f:
        try:
            while True:
                line = next(f)
                if "User input:" in line:
                    n_atoms, at_sym, method, basis = parse_sp_input(f)
                elif "Standard Nuclear Orientation" in line:
                    atom_nums, coord = parse_atoms_coords(f, n_atoms)
                    if atom_nums is None:
                        return None
                    atomic_nums.append(atom_nums)
                    coords.append(coord)
                # energy
                elif "Convergence failure" in line:
                    return None
                elif "Cycle" in line and "Energy" in line:
                    energies.append(parse_energy(f))
                # forces
                elif "Gradient of SCF Energy" in line:
                    force = parse_forces(f, n_atoms)
                    if force is not None:
                        forces.append(force)
                    else:
                        print(filename, " contains unparsed forces")
                # dipoles
                elif "Dipole Moment (Debye)" in line:
                    dipoles.append(parse_dipole(f))
                # quadrupoles
                elif "Quadrupole Moments (Debye-Ang)" in line:
                    quadrupoles.append(parse_quadrupole(f))
                # charges
                elif "Ground-State Mulliken Net Atomic Charges" in line:
                    charges.append(parse_charges(f, n_atoms))
                elif "Optimization Cycle" in line:
                    if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                            quadrupoles) == len(charges):
                        if line.split()[-1] == "1":
                            try:
                                atoms = [Atom(at_num) for at_num in atomic_nums[0]]
                                mset = MoleculeSet(atoms)
                                opt_trajectory = []
                            except:
                                print("Error making MSet for ", filename)
                                return None
                        if len(energies) > len(opt_trajectory):
                            opt_trajectory.append(build_new_geom(atomic_nums, coords, energies, forces, dipoles,
                                                                 quadrupoles, charges,
                                                                 method, basis))
                elif "OPTIMIZATION CONVERGED" in line:
                    mset.trajectories[".".join((method, basis, "opt"))] = opt_trajectory
                    return mset
        except StopIteration:
            print("Hit EOF on ", filename)
            try:
                mset.trajectories[".".join((method, basis, "opt"))] = opt_trajectory
                return mset
            except UnboundLocalError:
                print("No MSet built for ", filename)
                return None


def collect_meta_less50_data():
    for i in range(1, 10):
        for mol in glob.glob("/home/jeherr/tensorchem/data/chemspider_data/outputs/meta/less50/" + str(i) + "/*"):
            if os.path.isdir(mol):
                outfiles = glob.glob(os.path.join(mol, "*.out"))
                if len(outfiles) == 0:
                    print("No output files found in ", mol)
                    continue
                mol_mset = read_multi_sp_data(outfiles)
                head, name = os.path.split(mol)
                if mol_mset is not None:
                    mol_mset.save(filename=os.path.join("/mnt/sdb1/jeherr/chemspider_data/expanded_msets/meta/less_50/",
                                                        ".".join((name, "mset"))))


def collect_meta_data():
    for mol in glob.glob("/home/jeherr/tensorchem/data/chemspider_data/outputs/meta/crc/*.out") + glob.glob(
            "/home/jeherr/tensorchem/data/chemspider_data/outputs/meta/gigan/*.out") + glob.glob(
            "/home/jeherr/tensorchem/data/chemspider_data/outputs/meta/zerg/*.out"):
        head, name = os.path.split(mol)
        mol_mset = read_sp_data(mol)
        if mol_mset is not None:
            mol_mset.save(filename=os.path.join("/mnt/sdb1/jeherr/chemspider_data/expanded_msets/meta/first_runs/",
                                                ".".join((name[:-4], "mset"))))


def collect_aimd_data():
    for i in range(1, 5):
        for mol in glob.glob(f"/home/jeherr/tensorchem/data/chemspider_data/outputs/aimd/{i}/*.out"):
            head, name = os.path.split(mol)
            mol_mset = read_aimd_data(mol)
            if mol_mset is not None:
                mol_mset.save(filename=os.path.join("/mnt/sdb1/jeherr/chemspider_data/expanded_msets/aimd/",
                                                    ".".join((name[:-4], "mset"))))


def collect_opt_data():
    for i in range(5, 10):
        for mol in glob.glob(f"/home/jeherr/tensorchem/data/chemspider_data/outputs/opt/{i}/*.out"):
            head, name = os.path.split(mol)
            mol_mset = read_opt_data(mol)
            if mol_mset is not None:
                mol_mset.save(filename=os.path.join("/mnt/sdb1/jeherr/chemspider_data/expanded_msets/opt/",
                                                    ".".join((name[:-4], "mset"))))


#collect_aimd_data()
#collect_meta_less50_data()
#collect_meta_data()
#collect_opt_data()
