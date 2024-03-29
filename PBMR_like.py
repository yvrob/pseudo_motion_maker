from copy import deepcopy


## Simulation
neutron_pop = [1000000, 100, 5]
power_Wth = 400e6
acelib_path = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7u.xsdata"
dec_lib_path = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7.dec"
nfy_lib_path = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7.nfy"
opti = 1
ures = 1
additional_input_lines = ["set dd 2", "set gcu -1", "set pcc 0", "set bumode 2 16", "set memfrac 0.99", "set repro 0", "set shbuf 0 0", "set outp 10000000", ]
mpitasks = 20
parallel_creation = True
dep_output = "../one_peb9.8/input_one_pebble_dep.m"

# Core
Zmin = 0
Rin = 200 / 2  # inner radius, 0 if no inner radius
Rout = 370 / 2  # 120 # core radius
H = 1100  # 310 # core height
refl_thickness = 90
shape = "cyl"  # can be cyl or whatever (will be cuboid if not cyl)
quality = 6000  # pixels for plots

## Pebble bed
target_PF = 0.6  # packing fraction to reach
lattice_type = "fcc"  # can be hcp (not usable), fcc, sc

## Cycle
residence_time = 300
npasses = 10
average_discharge_burnup = 80
average_discharge_concentration = average_discharge_burnup*1.53e-06 # factor from PHYSOR2022 paper
direction = -1  # +1=up, -1=down
max_bu_step = 5  # days
ncycles_to_simulate = 5
threshold_type = "Cs137"

## Pebble
radii_pebble = [0, 2.5, 3.00]
# Parameters for iterative search with target packing fraction within volume
dist_pebbles_ini = radii_pebble[-1] * 5  # initial distance, huge on purpose
dist_pebbles_ini = 3.1387
dist_triso_ini = radii_pebble[-2]
mult_a = 1 - 1e-3  # initial multiplier for dist_pebbles, everytime the PF is too small
eps = 1e-3  # desired precision

## Triso
kernel_radius = 500e-4 / 2
buffer_radius = 93e-4 + kernel_radius
iPyC_radius = 38e-4 + buffer_radius
SiC_radius = 35e-4 + iPyC_radius
oPyC_radius = 40e-4 + SiC_radius
radii_triso = [
    kernel_radius,
    buffer_radius,
    iPyC_radius,
    SiC_radius,
    oPyC_radius,
]  # List of radii for the triso particles (cm) [fuel, buffer, inner PyC, SiC, outer PyC]
triso_PF = 15000

## Fuel
density = 10.4
density_a_or_m = "m"
enrich = 9.8e-2
n_U = 1
enrich_a_or_m = "m"
fuel_others_ZA_list = ["O16"]
fuel_others_n = [2]
fuel_temp = 900 + 300 #273.15

## Other materials
precreated_materials_path = None
# precreated_materials_path = "./materials"
reflector_material = "graphite"
inner_region_material = reflector_material
fuel_material = "fuel"
buffer_material = "buffer"
iPyC_material = "iPyC"
SiC_material = "SiC"
oPyC_material = "oPyC"
matrix_material = "matrix"
coolant_material = "helium"
central_graph_material = "central_graph"
graph_shell_material = "shell"
# to_check = [
#     reflector_material,
#     buffer_material,
#     iPyC_material,
#     SiC_material,
#     oPyC_material,
#     matrix_material,
#     coolant_material,
#     central_graph_material,
#     graph_shell_material,
# ]

other_materials_temp = 600  # K
helium_pressure = 60  # bar
helium_density = (
    48.14
    * helium_pressure
    / other_materials_temp
    * (1 + 0.4446 * helium_pressure / (other_materials_temp**1.2)) ** -1
) * 1e-3
helium = {
    "material_name": coolant_material,
    "ZA_list": ["He4"],
    "fractions_list": [1],
    "fractions_a_or_m": "a",
    "density": helium_density,
    "density_a_or_m": "m",
    "temperatureK": other_materials_temp,
    "moderator_name": None,
    "moderator_ZA": None,
}

reflector = {
    "material_name": reflector_material,
    "ZA_list": ["Cnat"],
    "fractions_list": [1],
    "fractions_a_or_m": "a",
    "density": 1.80,
    "density_a_or_m": "m",
    "temperatureK": other_materials_temp,
    "moderator_name": "others_therm",
    "moderator_ZA": "Cnat",
}


graph_shell = {
    "material_name": graph_shell_material,
    "ZA_list": ["Cnat"],
    "fractions_list": [1],
    "fractions_a_or_m": "a",
    "density": 1.75,
    "density_a_or_m": "m",
    "temperatureK": fuel_temp,
    "moderator_name": "triso_therm",
    "moderator_ZA": "Cnat",
}

central_graph = deepcopy(graph_shell)
central_graph["material_name"] = central_graph_material

matrix = deepcopy(graph_shell)
matrix["material_name"] = matrix_material

buffer = deepcopy(graph_shell)
buffer["material_name"] = buffer_material
buffer["density"] = 1.04

iPyC = deepcopy(buffer)
iPyC["material_name"] = iPyC_material
iPyC["density"] = 1.88

oPyC = deepcopy(buffer)
oPyC["material_name"] = oPyC_material
oPyC["density"] = 1.88

SiC = {
    "material_name": SiC_material,
    "ZA_list": ["Cnat", "Si28"],
    "fractions_list": [1, 1],
    "fractions_a_or_m": "a",
    "density": 3.15,
    "density_a_or_m": "m",
    "temperatureK": fuel_temp,
    "moderator_name": "triso_therm",
    "moderator_ZA": "Cnat",
}

other_materials = [
    helium,
    reflector,
    graph_shell,
    central_graph,
    matrix,
    buffer,
    iPyC,
    oPyC,
    SiC,
]
colors_materials = [
    [242, 242, 242],
    [253, 218, 236],
    [229, 216, 189],
    [255, 255, 204],
    [254, 217, 166],
    [222, 203, 228],
    [204, 235, 197],
    [179, 205, 227],
    [251, 180, 174],
]
for i, mat in enumerate(other_materials):
    mat["color"] = colors_materials[i]
fuel_color = [255, 0, 0]

## Thermal scattering libraries
triso_therm = {"name": "triso_therm", "temp": fuel_temp}
others_therm = {"name": "others_therm", "temp": other_materials_temp}
therm_scat_lib = [triso_therm, others_therm]


