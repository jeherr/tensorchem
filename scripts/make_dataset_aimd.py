import glob
import os
from ase.data import atomic_numbers
from tensorchem.dataset.molecule import MoleculeSet, Atom


def read_aimd_data(filename):
    n_atoms = None
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
                                    at_num = atomic_numbers[at_sym]
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
                # These next sections **SHOULD** only grab the first instance of these properties (i.e. before the AIMD
                # section starts)
                elif "Standard Nuclear Orientation" in line:
                    next(f)
                    next(f)
                    atom_nums, coord = [], []
                    for _ in range(n_atoms):
                        split_line = next(f).split()
                        at_sym = split_line[1]
                        at_num = atomic_numbers[at_sym]
                        atom_nums.append(at_num)
                        coord.append((float(split_line[2]), float(split_line[3]), float(split_line[4])))
                    atomic_nums.append(atom_nums)
                    coords.append(coord)
                # energy
                elif "Cycle" in line and "Energy" in line:
                    while True:
                        line = next(f)
                        if "Convergence criterion met" in line:
                            break
                    energies.append(float(line.split()[1]))
                # forces
                elif "Gradient of SCF Energy" in line:
                    force = []
                    while len(force) < n_atoms:
                        line = next(f)
                        new_atoms = len(line.split())
                        new_forces = [[] for _ in range(new_atoms)]
                        for _ in range(3):
                            line = next(f).split()
                            for atm, frc in enumerate(line[1:]):
                                new_forces[atm].append(float(frc) / -0.529177208590000)
                        force += new_forces
                    forces.append(force)
                # dipoles
                elif "Dipole Moment (Debye)" in line:
                    line = next(f).split()
                    dipole = [float(comp) for comp in line[1::2]]
                    dipoles.append(dipole)
                # quadrupoles
                elif "Quadrupole Moments (Debye-Ang)" in line:
                    quadrupole = []
                    for _ in range(2):
                        line = next(f).split()
                        quadrupole += [float(comp) for comp in line[1::2]]
                    quadrupoles.append(quadrupole)
                # charges
                elif "Ground-State Mulliken Net Atomic Charges" in line:
                    charge = []
                    for _ in range(3):
                        next(f)
                    for _ in range(n_atoms):
                        line = next(f).split()
                        charge.append(float(line[-1]))
                    charges.append(charge)
                elif "AB INITIO MOLECULAR DYNAMICS" in line:
                    atoms = [Atom(at_num) for at_num in atomic_nums[0]]
                    mset = MoleculeSet(atoms)
                    aimd_trajectory = []
                    mol_labels = {".".join((method, basis, "potential")): energies[0],
                                  ".".join((method, basis, "dipole")): dipoles[0],
                                  ".".join((method, basis, "quadrupole")): quadrupoles[0]}
                    atom_labels = {".".join((method, basis, "gradients")): forces[0],
                                   ".".join((method, basis, "charges")): charges[0]}
                    aimd_trajectory.append(mset.build_geom(coords[0], mol_labels, atom_labels))
                    while True:
                        line = next(f)
                        if "TIME STEPS COMPLETED" in line:
                            if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                                    quadrupoles) == len(charges):
                                mol_labels = {".".join((method, basis, "potential")): energies[-1],
                                              ".".join((method, basis, "dipole")): dipoles[-1],
                                              ".".join((method, basis, "quadrupole")): quadrupoles[-1]}
                                atom_labels = {".".join((method, basis, "gradients")): forces[-1],
                                               ".".join((method, basis, "charges")): charges[-1]}
                                aimd_trajectory.append(mset.build_geom(coords[-1], mol_labels, atom_labels))
                                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                                return mset
                            else:
                                print("Error reading AIMD trajectory at last step", filename, "Returning MSet before this step")
                                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                                return mset
                        elif "TIME STEP" in line:
                            time_step = int(line.split()[2].lstrip("#"))
                            time_au = float(line.split()[5])
                            time_fs = float(line.split()[8])
                            if len(atomic_nums) == len(energies) == len(forces) == len(coords) == len(dipoles) == len(
                                    quadrupoles) == len(charges):
                                mol_labels = {".".join((method, basis, "potential")): energies[-1],
                                              ".".join((method, basis, "dipole")): dipoles[-1],
                                              ".".join((method, basis, "quadrupole")): quadrupoles[-1]}
                                atom_labels = {".".join((method, basis, "gradients")): forces[-1],
                                               ".".join((method, basis, "charges")): charges[-1]}
                                aimd_trajectory.append(mset.build_geom(coords[-1], mol_labels, atom_labels))
                            else:
                                print("Error reading AIMD trajectory at step", time_step, filename, "Returning MSet before "
                                                                                                    "this step")
                                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                                return mset
                        elif "Standard Nuclear Orientation" in line:
                            next(f)
                            next(f)
                            atom_nums, coord = [], []
                            for _ in range(n_atoms):
                                split_line = next(f).split()
                                at_sym = split_line[1]
                                at_num = atomic_numbers[at_sym]
                                atom_nums.append(at_num)
                                coord.append((float(split_line[2]), float(split_line[3]), float(split_line[4])))
                            atomic_nums.append(atom_nums)
                            coords.append(coord)
                        # energy
                        elif "Cycle" in line and "Energy" in line:
                            while True:
                                line = next(f)
                                if "Convergence criterion met" in line:
                                    break
                            energies.append(float(line.split()[1]))
                        # forces
                        elif "Gradient of SCF Energy" in line:
                            force = []
                            while len(force) < n_atoms:
                                line = next(f)
                                new_atoms = len(line.split())
                                new_forces = [[] for _ in range(new_atoms)]
                                for _ in range(3):
                                    line = next(f).split()
                                    for atm, frc in enumerate(line[1:]):
                                        new_forces[atm].append(float(frc) / -0.529177208590000)
                                force += new_forces
                            forces.append(force)
                        # dipoles
                        elif "Dipole Moment (Debye)" in line:
                            line = next(f).split()
                            dipole = [float(comp) for comp in line[1::2]]
                            dipoles.append(dipole)
                        # quadrupoles
                        elif "Quadrupole Moments (Debye-Ang)" in line:
                            quadrupole = []
                            for _ in range(2):
                                line = next(f).split()
                                quadrupole += [float(comp) for comp in line[1::2]]
                            quadrupoles.append(quadrupole)
                        # charges
                        elif "Ground-State Mulliken Net Atomic Charges" in line:
                            charge = []
                            for _ in range(3):
                                next(f)
                            for _ in range(n_atoms):
                                line = next(f).split()
                                charge.append(float(line[-1]))
                            charges.append(charge)
        except StopIteration:
            print("Hit EOF on ", filename)
            try:
                mset.trajectories[".".join((method, basis, "aimd"))] = aimd_trajectory
                return mset
            except UnboundLocalError:
                print("No MSet built for ", filename)
                return None



for mol in glob.glob("/home/jeherr/tensorchem/data/chemspider_data/outputs/4/*.out"):
    head, name = os.path.split(mol)
    mol_mset = read_aimd_data(mol)
    if mol_mset is not None:
        mol_mset.save(filename=os.path.join("/mnt/sdb1/jeherr/chemspider_data/expanded_msets/aimd/4/",
                                        ".".join((name[:-4], "mset"))))
