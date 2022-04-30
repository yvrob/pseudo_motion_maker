import pandas as pd
import math
import numpy as np
from mixtures import *


def create_Serpent_fuel_material(
    enrich,
    n_U,
    others_ZA_list,
    others_n,
    enrich_a_or_m,
    density,
    density_a_or_m,
    temperatureK,
    total_volume=None,
    fuel_material="fuel",
    print_a_or_m="m",
    color=None,
    lib="endf",
):
    fuel = U_molecule(
        enrich,
        others_ZA_list,
        others_n,
        enrich_a_or_m,
        n_U,
        density,
        density_a_or_m,
        separate_U=True,
    )
    string = to_Serpent(
        fuel_material,
        fuel[0],
        temperatureK,
        total_volume,
        print_a_or_m=print_a_or_m,
        color=color,
        lib=lib,
    )
    return string


def create_Serpent_material(
    material_name,
    ZA_list,
    fractions_list,
    fractions_a_or_m,
    density,
    density_a_or_m,
    temperatureK,
    total_volume,
    moderator_name=None,
    moderator_ZA=None,
    print_a_or_m="m",
    color=None,
    lib="endf",
):
    if isinstance(ZA_list[0], str):
        names = ZA_list
        ZA_list = [name_to_ZA(i, separated=True) for i in ZA_list]
    else:
        names = [ZA_to_name(i) for i in ZA_list]
    atomic_masses = [atomic_mass(i[0], i[1]) for i in ZA_list]
    mat = mixture(
        fractions_list,
        atomic_masses,
        fractions_a_or_m=fractions_a_or_m,
        names=names,
        density=density,
        density_a_or_m="m",
    )
    string = to_Serpent(
        material_name,
        mat,
        temperatureK,
        total_volume,
        print_a_or_m=print_a_or_m,
        color=color,
        lib=lib,
        moderator_name=moderator_name,
        moderator_ZA=moderator_ZA,
    )

    return string


def create_graphite_therm(name, temp, lib="endf", graphite_lib="gre7"):
    if lib == "endf" and graphite_lib == "gre7":
        tmp_list = np.array([294, 400, 500, 600, 700, 800, 1000, 1200, 1600, 1999])
        markers_list = np.array(
            ["00", "04", "08", "12", "16", "18", "20", "22", "24", "26"]
        )
    idx = tmp_list[tmp_list <= temp].argmax()
    if temp == tmp_list[idx]:
        string = "therm\t{}\t{}.{}t\n".format(name, graphite_lib, markers_list[idx])
    else:
        string = "therm\t{}\t{}\t{}.{}t\t{}.{}t\n".format(
            name,
            temp,
            graphite_lib,
            markers_list[idx],
            graphite_lib,
            markers_list[idx + 1],
        )
    return string


def write_graphite_therm(
    path, name, temp, lib="endf", own_file=False, graphite_lib="gre7"
):
    string = create_graphite_therm(name, temp, lib=lib, graphite_lib=graphite_lib)
    if own_file:
        with open(path, "w") as f:
            f.write(string)
    else:
        with open(path, "a") as f:
            f.write(string + "\n")


def write_fuel_Serpent(
    path,
    enrich,
    n_U,
    others_ZA_list,
    others_n,
    enrich_a_or_m,
    density,
    density_a_or_m,
    temperatureK,
    total_volume=None,
    fuel_material="fuel",
    print_a_or_m="m",
    color=None,
    lib="endf",
    own_file=False,
):

    string = create_Serpent_fuel_material(
        enrich,
        n_U,
        others_ZA_list,
        others_n,
        enrich_a_or_m,
        density,
        density_a_or_m,
        temperatureK,
        total_volume,
        fuel_material=fuel_material,
        print_a_or_m=print_a_or_m,
        color=color,
        lib=lib,
    )
    if own_file:
        with open(path, "w") as f:
            f.write(string)
    else:
        with open(path, "a") as f:
            f.write(string + "\n")
    return True


def write_mat_Serpent(
    path,
    material_name,
    ZA_list,
    fractions_list,
    fractions_a_or_m,
    density,
    density_a_or_m,
    temperatureK,
    total_volume=None,
    print_a_or_m="m",
    moderator_name=None,
    moderator_ZA=None,
    color=None,
    lib="endf",
    own_file=False,
):

    string = create_Serpent_material(
        material_name,
        ZA_list,
        fractions_list,
        fractions_a_or_m,
        density,
        density_a_or_m,
        temperatureK,
        total_volume,
        moderator_name=moderator_name,
        moderator_ZA=moderator_ZA,
        print_a_or_m=print_a_or_m,
        color=color,
        lib=lib,
    )
    if own_file:
        with open(path, "w") as f:
            f.write(string)
    else:
        with open(path, "a") as f:
            f.write(string + "\n")

    return True


if __name__ == "__main__":
    # mass_dens = 10.6
    # enrich = 9.5e-2
    # wt_enrich = True
    # n_U = 1
    # n_O = 2
    # n_C = 0
    # verbose = True
    # fuel_temp = 1173
    # Npebbles = 250000
    # r_pebble = 2
    # fuel_volume = Npebbles * 4 / 3 * r_pebble**3
    # string = create_Serpent_fuel_material(
    #     mass_dens,
    #     enrich,
    #     n_U,
    #     n_O,
    #     n_C,
    #     fuel_temp,
    #     fuel_volume,
    #     verbose=True,
    #     wt_enrich=True,
    # )
    # print(string)

    fuel_densities = mixture([0.2, 0.8], [235, 238], "w")
