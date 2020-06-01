import glob
import os
import json
from ase.data import atomic_numbers

for mol in glob.glob("/Users/johnherr/tensorchem/data/chemspider_chno_opt/*.out"):
    with open(mol, "r") as f:
        lines = f.readlines()
    head, mol = os.path.split(mol)
    atoms = None

    for line in lines:
        if atoms is not None:
            if "$end" in line:
                break
            else:
                atoms += 1
        if "$molecule" in line:
            atoms = -1

    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []
    converged = False

    try:
        for i, line in enumerate(lines):
            if "**  OPTIMIZATION CONVERGED  **" in line:
                converged = True
                break
            # atomic numbers
            if "Standard Nuclear Orientation" in line:
                atom_nums = []
                for j in range(atoms - 1):
                    atom_nums.append(atomic_numbers[lines[i + 3 + j].split()[1]])
                atomic_nums.append(atom_nums)
            # coordinates
            if "Standard Nuclear Orientation" in line:
                xyz = []
                for j in range(atoms - 1):
                    xyz.append([float(lines[i + 3 + j].split()[2]), float(lines[i + 3 + j].split()[3]),
                                float(lines[i + 3 + j].split()[4])])
                coords.append(xyz)
            # energies
            if "Convergence criterion met" in line:
                energies.append(float(line.split()[1]))
            # forces
            if "Gradient of SCF Energy" in line:
                k = 0
                l = 0
                force = []
                for j in range(1, atoms):
                    force.append([float(lines[i + k + 2].split()[l + 1]) / -0.529177208590000,
                                  float(lines[i + k + 3].split()[l + 1]) / -0.529177208590000,
                                  float(lines[i + k + 4].split()[l + 1]) / -0.529177208590000])
                    l += 1
                    if (j % 6) == 0:
                        k += 4
                        l = 0
                forces.append(force)
            # dipoles
            if "Dipole Moment (Debye)" in line:
                dipole = []
                for j in range(3):
                    dipole.append(float(lines[i + 1].split()[1 + j * 2]))
                dipoles.append(dipole)
            # quadrupoles
            if "Quadrupole Moments (Debye-Ang)" in line:
                quadrupole = []
                for j in range(3):
                    quadrupole.append(
                        [float(lines[i + 1].split()[1 + j * 2]), float(lines[i + 2].split()[1 + j * 2])])
                quadrupoles.append(quadrupole)
            # charges
            if "Ground-State Mulliken Net Atomic Charges" in line:
                charge = []
                for j in range(atoms - 1):
                    charge.append(float(lines[i + j + 4].split()[2]))
                charges.append(charge)
            #print(len(atomic_nums), len(energies), len(forces), len(dipoles), len(quadrupoles), len(charges))
    except Exception as Ex:
        print(Ex, mol)
        continue

    if converged:
        if len(atomic_nums) == len(coords) == len(energies) == len(dipoles) == len(quadrupoles) == len(charges) > 0:
            for i, atomic_num_1 in enumerate(atomic_nums):
                for atomic_num_2 in atomic_nums[i+1:]:
                    if tuple(atomic_num_1) != tuple(atomic_num_2):
                        print("Atomic numbers are not the same between frames")
                        print(mol)
                        exit(0)
            mset_data = {"atomic_numbers": atomic_nums[len(atomic_nums) - 1]}
            geometries = [{"coordinates": coord} for coord in coords]
            for i, geom in enumerate(geometries):
                geom.update({"labels": {"wb97x-d.6-311gss.energy": energies[i],
                             "wb97x-d.6-311gss.forces": forces[i],
                             "wb97x-d.6-311gss.dipole": dipoles[i],
                             "wb97x-d.6-311gss.quadrupole": quadrupoles[i],
                             "wb97x-d.6-311gss.mulliken_charges": charges[i]}})
            mset_data.update({"geometries": geometries,
                              "identifiers": {"chemspider_id": mol[:-4]}})
            with open(os.path.join("/Users/johnherr/tensorchem/data/", ".".join((mol[:-4], "mset"))), "w") as f2:
                json.dump(mset_data, f2)
        else:
            print(mol)
            print(len(atomic_nums), len(energies), len(forces), len(dipoles), len(quadrupoles), len(charges))
            exit(0)
    else:
        print("Optimization did not converge for "+mol)
