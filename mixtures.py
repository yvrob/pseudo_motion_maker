"""
Functions for calculating mixtures and molecules

MAIN FUNCTIONS
- mixture: creates table with quantities for a given mixture of elements
mixture(
    fractions, <- list of quantities, does not have to sum up to one
    atomic_masses, <- list of atomic masses associated with quantities
    fractions_a_or_m, <- "m" if quantities are massic, "a" if quantities are atomic
    names=[], <- only for displaying, can be empty, or names linked to the quantities
    density, [optional] <- only if given, made to calculate not only wt/at fractions but also quantities or densities
    density_a_or_m <- "m" if density is massic, "a" if density is atomic
)
Note: masses in g. If density per cm^3 -> calculates in g/cm^3 and at/cm^3, yf density per m^3 -> calculates in g/m^3 and at/m^3, etc.
Example: 20 at% enriched uranium with mass density of 10.6 g/cm^3: 
    mix = mixture([0.2, 0.8], [atomic_mass(92,235), atomic_mass(92,238)], "a", ["U235", "U238"], 10.6, "m")

- molecule: creates a table with quantities for a molecule
molecule(
    ZA_list, <- list of [protons, nucleons] in molecule, natural if nucleons=0
    n_list, <- list of stoechiometric numbers for molecule
    density, [optional] <- only if given, made to calculate not only wt/at fractions but also quantities or densities
    density_a_or_m <- "m" if density is massic, "a" if density is atomic
)
Example: CO2 with mass density of 1.562 g/cm^3: 
    mol = molecule([[6,0],[8,0]], [1, 2], 1.562, "m")

- U_molecule: creates a table with quantities for a molecule containing enriched uranium
U_molecule(
    enrich, <- uranium enrichment (fraction)
    others_ZA_list, <- list of [protons, nucleons] in molecule except U, natural if nucleons=0
    others_n, <- list of stoechiometric numbers for molecule except U
    enrich_a_or_m <-"m" if enrichment is massic, "a" if enrichment is atomic
    n_U=1, <- stoechiometric number for uranium
    density, [optional] <- only if given, made to calculate not only wt/at fractions but also quantities or densities
    density_a_or_m, <- "m" if density is massic, "a" if density is atomic
    separate_U <- separate U235 and U238 in molecule table or keep U
)
Note: returns both the information about the molecule and uranium itself

Example: UO2 with mass density of 10.8 g/cm^3 and 20 wt% enrichment:
    Umol, U = U_molecule(0.2, [[8, 0]], [2], "m", 1, 10.8, "m")

OTHER FUNCTIONS
- Z_to_element: from number of protons, gives the name of the element
- element_to_Z: from element name (U, H, O, ...), gives number of protons
- atomic_mass: finds in NIST database the atomic (molar) mass from Z and A, if A<=0, finds natural atomic mass
- m_to_N: calculates from the mass in g and molar mass the atomic quantities in at (can be densities)
- N_to_m: calculates from the atomic quantities in at and molar mass the mass in g (can be densities)
- mixture_molar_mass: shorter way to extract molar mass of a mixture with fraction and atomic masses
"""

import numpy as np
import pandas as pd
import urllib.request
import math

# fmt: off
def ZA_to_name(ZA):
    ZA = str(ZA)
    # Case where ZA (+ temperature card) -> name
    ZA = int(ZA.split('.')[0])
    Z = int(ZA/1000)
    A = int(ZA-Z*1000)
    if A>=300:
        suffix = 'm'
        if Z > 80:
            A -= 100
        else:
            A -= 200
    else:
        suffix = ''
    if A == 0:
        A = 'nat'
    element = Z_to_element(Z)
    name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
    return name

def ZAI_to_name(ZAI):
    # Case where ZAI is given -> name
    if str(ZAI)[-1] == '0':
        ZA = int(ZAI)/10
        suffix = ''
    else:
        ZA = (int(ZAI) - int(str(ZAI)[-1]))/10
        suffix = 'm'

    Z = int(ZA/1000)
    if Z > 120:
        print('nan')
    element = Z_to_element(Z)
    A = int(ZA-Z*1000)
    if A == 0:
        A = 'nat'
    name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
    return name

def name_to_ZAI(name, separated=False):
    # Case where name is given -> ZAI
    if name[-1] == 'm':
        name = name[:-1]
        I = '1'
    else:
        I = '0'
    
    if name[-3:] == 'nat':
        element = name[:-3]
        Z = str(element_to_Z(element))
        if not separated:
            return int(Z+'000')
        else:
            return int(Z), 0

    i = len(name)-1
    while name[i:].isdigit():
        i-=1
    A = name[i+1:].zfill(3)
    element = name[:i+1]
    Z = str(element_to_Z(element))
    ZAI = int(Z+A+I)
    if not separated:
        return int(ZAI)
    else:
        return int(Z),int(A),int(I)
    
def name_to_ZA(name, separated=False):
    # Case where name is given -> ZA
    if name[-1] == 'm':
        name = name[:-1]
        isomer = True
    else:
        isomer = False
    if name[-3:] == 'nat':
        element = name[:-3]
        Z = str(element_to_Z(element))
        if not separated:
            return int(Z+'000')
        else:
            return int(Z), 0

    i = len(name)-1
    while name[i:].isdigit():
        i-=1
    A = name[i+1:].zfill(3)
    element = name[:i+1]
    Z = str(element_to_Z(element))
    if isomer:
        if int(A) < 100:
            A = str(int(A)+200)
        else:
            A = str(int(A)+100)
    ZA = Z+A
    if not separated:
        return int(ZA)
    else:
        return int(Z),int(A)

def Z_to_element(Z):
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
    return elements[Z-1]

def element_to_Z(element):
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
    return elements.index(element)+1
# fmt: on


def atomic_mass(Z, A):
    try:
        uf = urllib.request.urlopen(
            "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele={}&ascii=ascii2".format(
                Z
            )
        )
        if A <= 0:
            natural = True
        else:
            natural = False

        lines = uf.read().split(b"\n")
        lines = [line.decode("utf-8") for line in lines]

        for i, line in enumerate(lines):
            if natural and "Standard Atomic Weight" in line:
                extracted_text = line.split("Standard Atomic Weight = ")[1].split("(")[
                    0
                ]
                if extracted_text.isdigit():
                    M = float(extracted_text)
                else:
                    extracted_text = extracted_text.strip("][").split(",")
                    M = np.mean(np.array(extracted_text).astype(float))
                return M

            elif line == "Mass Number = {}".format(A):
                line_of_interest = lines[i + 1]
                extracted_text = line_of_interest.split("Relative Atomic Mass = ")[
                    1
                ].split("(")[0]
                M = float(extracted_text)
                return M
    except:
        return A
    raise Exception("Could not find isotope molar mass for Z={} and A={}".format(Z, A))


def m_to_N(m, atomic_mass):
    avogadro = 6.0221408e23
    n = m / atomic_mass
    N = n * avogadro
    return N


def N_to_m(N, atomic_mass):
    avogadro = 6.0221408e23
    n = N / avogadro
    m = n * atomic_mass
    return m


def mixture(
    fractions,
    atomic_masses,
    fractions_a_or_m="m",
    names=[],
    density=None,
    density_a_or_m="m",
):
    avogadro = 6.0221408e23

    try:
        test = fractions[0]
    except:
        fractions = [fractions]
        atomic_masses = [atomic_masses]
    if not isinstance(density, type(None)):
        if density_a_or_m == "w":
            density_a_or_m = "m"
        if density_a_or_m not in ["m", "a"]:
            raise Exception("Choose atomic (a) or mass (m) nature for density")

    if fractions_a_or_m == "w":
        fractions_a_or_m = "m"
    if fractions_a_or_m not in ["m", "a"]:
        raise Exception("Choose atomic (a) or mass (m) nature for fractions")

    atoms = []
    masses = []
    for i in range(len(fractions)):
        if fractions_a_or_m == "m":
            masses.append(fractions[i])
            atoms.append(m_to_N(fractions[i], atomic_masses[i]))
        elif fractions_a_or_m == "a":
            atoms.append(fractions[i])
            masses.append(N_to_m(fractions[i], atomic_masses[i]))
    mix = pd.DataFrame()
    mix["N"] = atoms
    mix["m"] = masses
    mix["M"] = atomic_masses
    if len(names) == len(fractions):
        mix["Element"] = names
    else:
        mix["Element"] = atomic_masses
    mix = mix.set_index("Element")

    mix["at_fraction"] = mix["N"] / mix["N"].sum()
    mix["at_percent"] = mix["at_fraction"] * 100
    mix["wt_fraction"] = mix["m"] / mix["m"].sum()
    mix["wt_percent"] = mix["wt_fraction"] * 100
    mix["stoech"] = np.divide(mix["at_fraction"], mix["at_fraction"].min())

    if not isinstance(density, type(None)):
        if density_a_or_m == "m":
            mix["m"] = density * mix["wt_fraction"]
            mix["N"] = mix["m"] / mix["M"] * avogadro
        elif density_a_or_m == "a":
            mix["N"] = density * mix["at_fraction"]
            mix["m"] = mix["N"] / avogadro * mix["M"]
    else:
        mix = mix.drop(columns=["m", "N"])

    M = 1 / np.sum((mix["wt_fraction"] / mix["M"]))
    mix.loc["Mixture"] = mix.sum(numeric_only=True, axis=0)
    mix.loc["Mixture"]["M"] = M
    return mix


def mixture_molar_mass(fractions, atomic_masses, fractions_a_or_m="m"):
    mix = mixture(fractions, atomic_masses, fractions_a_or_m)
    return mix["M"]["Mixture"]


def molecule(
    ZA_list,
    n_list,
    density=None,
    density_a_or_m="m",
):
    avogadro = 6.0221408e23

    if not isinstance(density, type(None)):
        if density_a_or_m == "w":
            density_a_or_m = "m"
        if density_a_or_m not in ["m", "a"]:
            raise Exception("Choose atomic (a) or mass (m) nature for density")

    names = []
    M = []
    n = []
    if not isinstance(ZA_list[0], (list, tuple, np.ndarray)):
        ZA_list = [name_to_ZA(i, separated=True) for i in ZA_list]
    for i in range(len(ZA_list)):
        Z = ZA_list[i][0]
        A = ZA_list[i][1]
        names.append("{}{}".format(Z_to_element(Z), A if A > 0 else "nat"))
        M.append(atomic_mass(Z, A))
        n.append(n_list[i])
    mol = mixture(n, M, "a", names, density=density, density_a_or_m=density_a_or_m)
    mol = mol.drop(index="Mixture")
    mol["stoech"] = n_list
    M_mol = np.sum(np.multiply(M, n))
    mol.loc["Molecule"] = mol.sum(numeric_only=True, axis=0)
    mol.loc["Molecule"]["M"] = M_mol
    return mol


def U_molecule(
    enrich,
    others_ZA_list,
    others_n,
    enrich_a_or_m="m",
    n_U=1,
    density=None,
    density_a_or_m="m",
    separate_U=True,
):
    avogadro = 6.0221408e23

    if enrich_a_or_m == "w":
        enrich_a_or_m = "m"
    if enrich_a_or_m not in ["m", "a"]:
        raise Exception("Choose atomic (a) or mass (m) nature for enrichment")
    U = mixture(
        [enrich, 1 - enrich],
        [atomic_mass(92, 235), atomic_mass(92, 238)],
        enrich_a_or_m,
        names=["U235", "U238"],
    )
    M_U = U["M"]["Mixture"]
    names = ["U"]
    M = [M_U]
    n = [n_U]
    if not isinstance(others_ZA_list[0], (list, tuple, np.ndarray)):
        others_ZA_list = [name_to_ZA(i, separated=True) for i in others_ZA_list]
    for i in range(len(others_ZA_list)):
        Z = others_ZA_list[i][0]
        A = others_ZA_list[i][1]
        names.append("{}{}".format(Z_to_element(Z), A if A > 0 else "nat"))
        M.append(atomic_mass(Z, A))
        n.append(others_n[i])
    Umol = mixture(n, M, "a", names, density=density, density_a_or_m=density_a_or_m)

    if separate_U:
        Umol = U.append(Umol, ignore_index=False, sort=False)
        Umol.loc["U235"]["N"] = Umol.loc["U"]["N"] * Umol.loc["U235"]["at_fraction"]
        Umol.loc["U235"]["at_fraction"] = (
            Umol.loc["U"]["at_fraction"] * Umol.loc["U235"]["at_fraction"]
        )
        Umol.loc["U235"]["m"] = Umol.loc["U"]["m"] * Umol.loc["U235"]["wt_fraction"]
        Umol.loc["U235"]["wt_fraction"] = (
            Umol.loc["U"]["wt_fraction"] * Umol.loc["U235"]["wt_fraction"]
        )

        Umol.loc["U238"]["N"] = Umol.loc["U"]["N"] * Umol.loc["U238"]["at_fraction"]
        Umol.loc["U238"]["at_fraction"] = (
            Umol.loc["U"]["at_fraction"] * Umol.loc["U238"]["at_fraction"]
        )
        Umol.loc["U238"]["m"] = Umol.loc["U"]["m"] * Umol.loc["U238"]["wt_fraction"]
        Umol.loc["U238"]["wt_fraction"] = (
            Umol.loc["U"]["wt_fraction"] * Umol.loc["U238"]["wt_fraction"]
        )
        Umol = Umol.drop(index="U")
    Umol = Umol.drop(index="Mixture")
    if separate_U:
        Umol["stoech"] = [
            U.loc["U235"]["at_fraction"],
            U.loc["U238"]["at_fraction"],
        ] + others_n
    else:
        Umol["stoech"] = n
    Umol.loc["Molecule"] = Umol.sum(numeric_only=True, axis=0)
    return Umol, U


def mix_two(table1, table2):
    mix = pd.DataFrame().append(table1.add(table2, fill_value=0), sort=False)
    if "Molecule" in mix.index:
        name = "Molecule"
    elif "Mixture" in mix.index:
        name = "Mixture"

    mix = mix.drop(index=name)
    mix["at_fraction"] /= mix["at_fraction"].sum()
    mix["wt_fraction"] /= mix["wt_fraction"].sum()
    mix["at_percent"] /= mix["at_percent"].sum()
    mix["wt_percent"] /= mix["wt_percent"].sum()
    for i in mix.index:
        try:
            mix.loc[i]["M"] = table1.loc[i]["M"]
            mix.loc[i]["stoech"] = table1.loc[i]["stoech"]
        except:
            mix.loc[i]["M"] = table2.loc[i]["M"]
            mix.loc[i]["stoech"] = table2.loc[i]["stoech"]
    M = np.sum(np.multiply(mix["M"], mix["stoech"]))
    mix.loc[name] = mix.sum(numeric_only=True, axis=0)
    mix.loc[name]["M"] = M
    idx = [i for i in mix.index if i != name]
    mix = pd.DataFrame().append(mix, sort=False).loc[idx + [name]]
    return mix


def to_Serpent(
    name,
    mix,
    temperatureK,
    volume=None,
    moderator_name=None,
    moderator_ZA=None,
    additional_lines=[],
    print_a_or_m="m",
    color=None,
    lib="endf",
):
    if print_a_or_m == "w":
        print_a_or_m = "m"
    if print_a_or_m not in ["m", "a"]:
        raise Exception("Choose atomic (a) or mass (m) nature for density")
    if lib == "endf":
        temp_marker = 3 * math.floor(temperatureK / 300)
        if temp_marker == temperatureK / 100:
            string_temp = ""
        else:
            string_temp = "\ttmp\t{}".format(temperatureK)
        temp_marker = "%02d" % temp_marker

    if isinstance(volume, type(None)):
        burn = False
    else:
        burn = True
    if isinstance(moderator_name, type(None)):
        mod = False
    else:
        mod = True
    if "Mixture" in mix.index:
        field = "Mixture"
    elif "Molecule" in mix.index:
        field = "Molecule"
    string = "mat\t{}\t{:.6E}{}".format(
        name,
        -mix.loc[field]["m"] if print_a_or_m == "m" else mix.loc[field]["n"],
        string_temp,
    )
    if isinstance(moderator_ZA, str):
        moderator_ZA = name_to_ZA(moderator_ZA)
    if mod:
        string += "\tmoder\t{}\t{}".format(moderator_name, moderator_ZA)
    if burn:
        string += "\tburn\t1\tvol\t{}".format(volume)
    if not isinstance(color, type(None)):
        string += "\trgb\t{} {} {}".format(color[0], color[1], color[2])
    string += "\t% {}dens = {:.6E}\n".format(
        "a" if print_a_or_m == "m" else "n",
        mix.loc[field]["N"] if print_a_or_m == "m" else mix.loc[field]["m"],
    )
    for i in mix.index:
        if i != field:
            string += "{}.{}c\t{:.6E}".format(
                name_to_ZA(i),
                temp_marker,
                -mix.loc[i]["m"] if print_a_or_m == "m" else mix.loc[i]["N"],
            )
            string += "\t% {}dens = {:.6E}\n".format(
                "a" if print_a_or_m == "m" else "n",
                mix.loc[i]["N"] if print_a_or_m == "m" else mix.loc[i]["m"],
            )
    for line in additional_lines:
        string += line + "\n"
    return string
