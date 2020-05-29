import glob, os
from ase.data import atomic_numbers

if not os.path.exists("/home/adriscoll/tensormol-jax/tensormol_jax/data/chemspider_data/"):
    os.makedirs("/home/adriscoll/tensormol-jax/tensormol_jax/data/chemspider_data/")
for mol in glob.glob("/home/adriscoll/tensormol-jax/tensormol_jax/data/chemspider_data/19021.out"):
    f=open(mol, "r")
    lines = f.readlines()
    head, mol = os.path.split(mol)
    atoms = None
	
    for line in lines:
        if (atoms != None):
            if line.count("$end") > 0:
                break
            else:
                atoms += 1
        if line.count("$molecule") > 0:
            atoms = -1

    atomic_nums = []
    coords = []
    energies = []
    forces = []
    dipoles = []
    quadrupoles = []
    charges = []

    try:	
        for i, line in enumerate(lines):
            if "**  OPTIMIZATION CONVERGED  **" in line:
                break
            #atomic numbers
            if "Standard Nuclear Orientation" in line:
                atom_nums = []
                for j in range(atoms-1):
                    atom_nums.append(atomic_numbers[lines[i+3+j].split()[1]])
                atomic_nums.append(atom_nums)
            #coordinates
            if "Standard Nuclear Orientation" in line:
                xyz = []
                for j in range(atoms-1):
                    xyz.append([lines[i+3+j].split()[2],lines[i+3+j].split()[3],lines[i+3+j].split()[4]])
                coords.append(xyz)
            #energies
            if "Convergence criterion met" in line:
                energies.append(line.split()[1])
            #forces
            if "Gradient of SCF Energy" in line:
                k = 0
                l = 0
                force = []
                for j in range(1, atoms):
                    force.append([str(float(lines[i+k+2].split()[l+1])/-0.529177208590000),
                            str(float(lines[i+k+3].split()[l+1])/-0.529177208590000),
                            str(float(lines[i+k+4].split()[l+1])/-0.529177208590000)])
                    l += 1
                    if (j % 6) == 0:
                        k += 4
                        l = 0
                    forces.append(force)
            #dipoles
            if "Dipole Moment (Debye)" in line:
                dipole = []
                for j in range(3):
                    dipole.append(str(float(lines[i+1].split()[1+j*2])))
                dipoles.append(dipole)
            #quadrupoles
            if "Quadrupole Moments (Debye-Ang)" in line:
                quadrupole = []
                for j in range(3):
                    quadrupole.append([str(float(lines[i+1].split()[1+j*2])),str(float(lines[i+2].split()[1+j*2]))])
                quadrupoles.append(quadrupole)
            #charges
            if line.count("Ground-State Mulliken Net Atomic Charges") > 0:
                charge = []
                for j in range(atoms-1):
                    charge.append(lines[i+j+4].split()[2])
                charges.append(charge)
    except Exception as Ex:
        print(Ex, mol)
        continue

    if len(atomic_nums) == len(coords) == len(energies) == len(dipoles) == len(quadrupoles) == len(charges) > 0:
        f2 = open("/home/adriscoll/tensormol-jax/tensormol_jax/data/chemspider_data/" + mol[:-4] + ".mset", "w")
        x = len(energies)-1
        y = len(coords[0])-1
        f2.write('{"atomic_number": ' + str(atomic_nums[len(atomic_nums)-1]) + ', ')
        f2.write('"coordinates": [')
        for i in range(x+1):
            f2.write('[')
            for j in range(y):
                f2.write("[" + str(coords[i][j][0]) + ", " + str(coords[i][j][1]) + ", " + str(coords[i][j][2]) + "], ")
            if i == x:
                f2.write('[' + str(coords[i][y][0]) + ', ' + str(coords[i][y][1]) + ', ' + str(coords[i][y][2]) + ']]], ')
            else:
                f2.write('[' + str(coords[i][y][0]) + ', ' + str(coords[i][y][1]) + ', ' + str(coords[i][y][2]) + ']], ')
        f2.write('"properties": {"energies": ' + str(float(energies[x])) + ', "forces": ' + str(forces[x]).replace("'", "") + ', "dipoles": '\
            + str(dipoles[x]).replace("'", "") + ', "quadrupoles": ' + str(quadrupoles[x]).replace("'", "") + ', "charges": ' + str(charges[x]).replace("'", "") + '}}')
        f2.close()
    else:
        print(mol)
