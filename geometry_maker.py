#%%
# %load_ext autoreload
# %autoreload 2

from packer import find_lattice
from pbed_creator import write_pbed
from Serpent_geometry_creator import write_geometry_file
from Serpent_materials_creator import (
    write_fuel_Serpent,
    write_mat_Serpent,
    write_graphite_therm,
)
from precreated_materials_importer import import_precreated
from Serpent_simulation_creator import write_input
from burnup_calculator import calc_initial_BU
from first_restart_creator import create_decomposed_restart_files, replicate_dd
from motion_sequence_creator import write_motion_sequence,  write_bustep, create_sequence_scheme
from plotter import *
from mixtures import name_to_ZAI

from glob import glob
import os
import shutil
from subprocess import call
import importlib
import warnings

#%% Input
# What to do
create_geometry    = True
create_Serpent     = True
calculate_BU_distr = True
create_restart     = True
calculate_motion   = True

## Paths and naming
path = "./PBMR_like"
name_pos = "fpb_pos"
universe_pebbles = "u_pebble"

filename_input = "PBMR_like"


# Secondary input
same_rows = True  # needed for discrete motion
partial_pebbles = False  # allow for pebbles crossing bounds in volume
plot = True

#%% Read input from file
globals().update(importlib.import_module(filename_input.replace(".py", "")).__dict__)

path_model = path + "/model"
os.makedirs(path_model, exist_ok=True)

if create_geometry:
    print('Creating geometry')
    #%% Iteratively find pebble lattice positions
    pebble_bed, dist_pebbles, PF, found = find_lattice(
        lattice_type,
        radii_pebble[-1],
        shape,
        Rout,
        H,
        target_PF,
        Rin,
        refl_thickness,
        dist_pebbles_ini,
        mult_a,
        eps,
        partial_pebbles=partial_pebbles,
        centeredR=True,
        same_rows=same_rows,
    )
    if not found:
        warnings.warn(
            "Did not find any good solution! Taking best with dist_pebbles={:.4f}, (PF={:.4f})".format(
                dist_pebbles * mult_a**2, PF
            )
        )
    pbed_written = write_pbed(
        path_model,
        name_pos,
        np.array(pebble_bed[["x", "y", "z"]]),
        np.array(pebble_bed["r_pebbles"]),
        universe_pebbles,
    )

if create_Serpent:
    print('Creating Serpent input files')
    #%% Write geometry to input files
    with open(path_model + "/materials", "w") as f:
        pass
    geometry_written = write_geometry_file(
        path_model,
        shape,
        Rin,
        Rout,
        H,
        Zmin,
        refl_thickness,
        radii_pebble,
        universe_pebbles,
        radii_triso,
        triso_PF,
        name_pos,
        inner_region_material,
        reflector_material,
        fuel_material,
        buffer_material,
        iPyC_material,
        SiC_material,
        oPyC_material,
        matrix_material,
        coolant_material,
        central_graph_material,
        graph_shell_material,
        plot=plot,
        quality=quality,
        dist_triso_ini=dist_triso_ini,
    )

    #%% Plot and extract numbers and volumes
    Npebbles = len(import_last(path_model, pattern='fpb_pos'))
    Ntrisos = len(import_last(path_model, pattern="trisos*.inp"))
    fuel_volume = (4 / 3 * np.pi * radii_triso[0] ** 3) * Ntrisos * Npebbles


    #%% Materials
    # Import precreated materials
    if precreated_materials_path != None:
        import_precreated(
            precreated_materials_path,
            path_model,
            check=to_check,
        )

    fuel_written = write_fuel_Serpent(
        path_model + "/materials",
        enrich,
        n_U,
        fuel_others_ZA_list,
        fuel_others_n,
        enrich_a_or_m,
        density,
        density_a_or_m,
        fuel_temp,
        fuel_volume,
        fuel_material=fuel_material,
        print_a_or_m="m",
        color=fuel_color,
        lib="endf",
        own_file=False,
    )

    if not isinstance(other_materials, type(None)):
        for mat in other_materials:
            mat_written = write_mat_Serpent(
                path_model + "/materials",
                material_name=mat["material_name"],
                ZA_list=mat["ZA_list"],
                fractions_list=mat["fractions_list"],
                fractions_a_or_m=mat["fractions_a_or_m"],
                density=mat["density"],
                density_a_or_m=mat["density_a_or_m"],
                temperatureK=mat["temperatureK"],
                print_a_or_m="m",
                moderator_name=mat["moderator_name"],
                moderator_ZA=mat["moderator_ZA"],
                color=mat["color"],
                lib="endf",
                own_file=False,
            )

    for therm in therm_scat_lib:
        write_graphite_therm(path_model + "/materials", therm["name"], therm["temp"])

    #%% Create simulation
    if threshold_type == 'burnup':
        threshold = average_discharge_burnup
    elif threshold_type == 'passes':
        threshold = npasses
    else:
        threshold_type = str(name_to_ZAI(threshold_type))
        threshold = average_discharge_concentration

    input_written = write_input(
        path_model,
        neutron_pop,
        power_Wth,
        acelib_path,
        dec_lib_path,
        nfy_lib_path,
        opti,
        ures,
        threshold_type,
        threshold,
        fuel_material,
        additional_input_lines,
        initial_restart=True,
        check=False,
    )
pebble_bed = import_last(path_model, pattern='fpb_pos')
if calculate_BU_distr:
    print('Creating initial BU distribution')

    #%% Calculate BU guess and apply comp = f(BU)
    pebble_bed = calc_initial_BU(
        pebble_bed, residence_time, npasses, average_discharge_burnup, H
    )


#%% Create first restart files
if create_restart:
    print('Creating initial restart file')
    from time import time
    t0 = time()
    pebble_bed = replicate_dd(pebble_bed, mpitasks)
    create_decomposed_restart_files(
        pebble_bed, mpitasks, parallel_creation, dep_output, path_model, name_out="first_compos"
    )
    print("{} s elapsed".format(time() - t0))

#%% Motion sequence
if calculate_motion:
    print('Calculating motion')
    sequence_scheme = create_sequence_scheme(pebble_bed, max_bu_step, residence_time, npasses, ncycles_to_simulate)
    write_bustep(path_model+'/input', pebble_bed, sequence_scheme, residence_time, npasses)
    write_motion_sequence(path_model, pebble_bed, direction, sequence_scheme, lattice_type, npasses, residence_time , randomize_reloading=True, plot=False)

