from glob import glob
import os
import numpy as np


def create_input(
    neutron_pop,
    power_Wth,
    acelib_path,
    dec_lib_path,
    nfy_lib_path,
    opti,
    ures,
    threshold_type,
    threshold,
    fuel_material="fuel",
    additional_lines=[],
    initial_restart=True,
):
    string = "include materials\n"
    string += "include geometry\n"

    string += "set seed 12345\n"
    string += "\nset memfrac 0.99\n"
    string += "\ndiv {} sep 4\n".format(fuel_material)
    string += "set his 1\n"
    string += 'set inventory "all"\n'
    if len(acelib_path) > 0:
        string += 'set acelib "{}"\n'.format(acelib_path)
    if len(dec_lib_path) > 0:
        string += 'set declib "{}"\n'.format(dec_lib_path)
    if len(nfy_lib_path) > 0:
        string += 'set nfylib "{}"\n'.format(nfy_lib_path)
    string += "set pop {}\n".format(" ".join(np.array(neutron_pop).astype(str)))
    string += "set power {}\n".format(power_Wth)
    string += "set opti {}\n".format(opti)
    string += "set ures {}\n".format(ures)
    string += "set printm 1 0.0\n"
    if initial_restart:
        string += 'set rfr idx 0 "./restart/first_compos.wrk"\n'

    string += 'set mixfile "./indices/input_mat_indices" {} threshold "{}" {} waste "waste_file"\n'.format(
        fuel_material, threshold_type, threshold
    )
    for i in additional_lines:
        string += i+'\n'

    return string


# def list_dependencies(path, input_file):
#     with open(path+'/'+input_file) as f:
#         lines = f.readlines()
#     files = []
#     for line in lines:
#         if line[: len("include")] == "include":
#             files.append(line.split()[1].replace('"', ""))
#         elif line[: len("pbed")] == "pbed":
#             files.append(line.split()[3].replace('"', ""))
#     if len(files) == 0:
#         return []
#     else:
#         for
#         return list_dependencies(path, file)

# def check_existence(input_file):
#     files = list_dependencies(input_file)
#     if len(files) == 0:
#         return files
#     for file in list_dependencies(input_file):
#         if
#         check_existence(input_file)
#     return os.path.exists(input_file)


def write_input(
    path,
    neutron_pop,
    power_Wth,
    acelib_path,
    dec_lib_path,
    nfy_lib_path,
    opti,
    ures,
    threshold_type,
    threshold,
    fuel_material="fuel",
    additional_lines=[],
    initial_restart=True,
    *,
    check=True
):
    with open(path + "/input", "w") as f:
        f.write(
            create_input(
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
                additional_lines,
                initial_restart,
            )
        )

    # if check:
    #     check_existence(path + "/input")
