#%%
from packer import find_lattice
from pbed_creator import write_pbed
import numpy as np
import os


def create_vessel(
    shape,
    Rin,
    Rout,
    H,
    Zmin=0,
    refl_thickness=0,
    graphite_material="graphite",
    inner_region_material="graphite",
):
    Zmin_core = Zmin + refl_thickness
    Zmax_core = Zmin_core + H
    Zmin_vessel = Zmin
    Zmax_vessel = Zmax_core + refl_thickness
    Rout_vessel = Rout + refl_thickness

    # Surfaces to create core vessel
    string = "%%---surf for core vessel\nsurf infinite inf\n"

    if shape == "cyl":
        string += "surf core_surf_out cylz 0 0 {} {} {}\n".format(
            Rout, Zmin_core, Zmax_core
        )
        string += "surf core_surf_in  cylz 0 0 {} {} {}\n".format(
            Rin, Zmin_core, Zmax_core
        )
    elif shape == "cube":
        string += "surf core_surf_out cuboid {} {} {} {} {} {}\n".format(
            -Rout, Rout, -Rout, Rout, Zmin_core, Zmax_core
        )
        string += "surf core_surf_in cuboid {} {} {} {} {} {}\n".format(
            -Rin, Rin, -Rin, Rin, Zmin_core, Zmax_core
        )

    if refl_thickness != 0:
        if shape == "cyl":
            string += """
%%---graphite reflector surfaces
surf vessel_surf_out cyl 0.0 0.0 {} {} {} % outer reflector surface
""".format(
                Rout_vessel, Zmin_vessel, Zmax_vessel
            )
        elif shape == "cube":
            string += """
%%---graphite reflector surfaces
surf vessel_surf_out cuboid {} {} {} {} {} {} % outer reflector surface 
""".format(
                -Rout_vessel,
                Rout_vessel,
                -Rout_vessel,
                Rout_vessel,
                Zmin_vessel,
                Zmax_vessel,
            )
        string += "%%---core vessel cells\n"
        if Rin != 0:
            string += "cell inner_region_cell 0 {} -core_surf_in\n".format(
                inner_region_material
            )
        string += """cell reflector_cell 0 {} core_surf_out -vessel_surf_out
cell out_cell 0 outside vessel_surf_out
""".format(
            graphite_material
        )
    else:
        string += """
%%--- core vessel cells
cell out_cell 0 outside vessel_surf_out
"""
    return string


def create_pbed_geometry(
    pos_file,
    radii_pebbles,
    Rin,
    universe_pebbles,
    coolant_material="Flibe",
    central_graph_material="central_graph",
    graph_shell_material="shell",
):
    string = "%%---surf for fuel pebbles\n"
    if radii_pebbles[0] != 0:
        string += "surf central_graphite sph 0 0 0 {} % internal graphite\n".format(
            radii_pebbles[0]
        )
    string += (
        "surf graphite_mat sph 0 0 0 {} % graphite matrix maximum radius\n".format(
            radii_pebbles[1]
        )
    )

    string += "\n%---cells for fuel pebbles\n"
    if radii_pebbles[0] != 0:
        string += "cell c_pebble_center {} {} -central_graphite\n".format(
            universe_pebbles, central_graph_material
        )

        string += "cell c_pebble_matrix {} fill trisos -graphite_mat central_graphite\n".format(
            universe_pebbles, central_graph_material
        )
    else:
        string += "cell c_pebble_matrix {} fill trisos -graphite_mat\n".format(
            universe_pebbles, central_graph_material
        )
    if radii_pebbles[2] != 0:
        string += "cell c_pebble_shell {} {} graphite_mat\n".format(
            universe_pebbles, graph_shell_material
        )

    # Create pebble ped and coolant universes
    string += "\n%%---pebble bed\n"
    string += 'pbed u_pb u_coolant "{}"\n'.format(pos_file)
    if Rin > 0:
        string += "cell c_pb 0 fill u_pb  -core_surf_out core_surf_in\n"
    else:
        string += "cell c_pb 0 fill u_pb  -core_surf_out\n"
    string += "cell c_coolant u_coolant {} -infinite\n\n".format(coolant_material)
    return string


def create_triso(
    radii_triso,
    radii_pebble,
    path_to_geometry,
    fuel_material="fuel",
    buffer_material="buffer",
    iPyC_material="iPyC",
    SiC_material="SiC",
    oPyC_material="oPyC",
    matrix_material="matrix",
    *,
    triso_file=None,
    triso_PF=None,
    dist_triso_ini = None
):
    if isinstance(dist_triso_ini, type(None)):
        dist_triso_ini = radii_pebble[1]

    if (isinstance(triso_PF, type(None)) and isinstance(triso_file, type(None))) or (
        not isinstance(triso_PF, type(None)) and not isinstance(triso_file, type(None))
    ):
        raise Exception("Can only import pbed file or (not and) lattice distance")

    string = "%---triso particles\n"
    string += "particle triso_particle\n{}\t{}\n".format(fuel_material, radii_triso[0])
    if len(radii_triso) == 5:
        string += "{}\t{}\n".format(buffer_material, radii_triso[1])
        string += "{}\t{}\n".format(iPyC_material, radii_triso[2])
        string += "{}\t{}\n".format(SiC_material, radii_triso[3])
        string += "{}\t{}\n".format(oPyC_material, radii_triso[4])
        string += "{}\n\n".format(matrix_material)

    if isinstance(triso_file, type(None)):
        triso_file = "trisos_{}.inp".format(triso_PF)
        df, dist, PF, found = find_lattice(
            "sc",
            radii_triso[-1],
            "sph",
            radii_pebble[1],
            0,
            triso_PF,
            radii_pebble[0],
            0,
            dist_triso_ini,
            1 - 1e-2,
            1e-3,
            partial_pebbles=False,
            centeredR=True,
            centeredZ=True,
            same_rows=False,
            verbose=False,
        )
        pbed_written = write_pbed(
            path_to_geometry,
            triso_file,
            np.array(df[["x", "y", "z"]]),
            np.array(df["r_pebbles"]),
            "triso_particle",
        )
    string += 'cell inf u_matrix {} -infinite\n'.format(matrix_material)
    string += 'pbed trisos u_matrix "{}"\n'.format(triso_file)
    return string


def create_plot(quality, shape, Rout, Rin, H, Zmin, radii_pebble):
    if shape == "cyl":
        delta_x = 2 * Rout
        delta_z = H
    elif shape == "sph":
        delta_x = 2 * Rout
        delta_z = 2 * Rout
    else:
        delta_x = Rout - Rin
        delta_z = H
    string = """%%---plotting geometry
plot 1 {0} {1}
plot 1 {0} {1} {2}
plot 1 {0} {1} {3}
plot 3 {1} {1} {4}
plot 3 {1} {1} {5}
plot 3 {1} {1} {6}\n\n""".format(
        int(quality * (delta_x / delta_z)),
        quality,
        radii_pebble[-1]/2,
        radii_pebble[-1],
        Zmin + H / 2,
        Zmin + H / 2 + radii_pebble[-1]/2,
        Zmin + H / 2 + radii_pebble[-1],
    )
    return string


def write_geometry_file(
    path,
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
    inner_region_material="graphite",
    reflector_material="graphite",
    fuel_material="fuel",
    buffer_material="buffer",
    iPyC_material="iPyC",
    SiC_material="SiC",
    oPyC_material="oPyC",
    matrix_material="matrix",
    coolant_material="Flibe",
    central_graph_material="central_graph",
    graph_shell_material="shell",
    *,
    plot=False,
    quality=1000,
    dist_triso_ini = None
):
    if plot:
        block_plot = create_plot(quality, shape, Rout, Rin, H, Zmin, radii_pebble)
    else:
        block_plot = ""
    block_geom = create_vessel(
        shape,
        Rin,
        Rout,
        H,
        Zmin,
        refl_thickness,
        reflector_material,
        inner_region_material,
    )
    block_pbed = create_pbed_geometry(
        name_pos,
        radii_pebble,
        Rin,
        universe_pebbles,
        coolant_material,
        central_graph_material,
        graph_shell_material,
    )
    block_triso = create_triso(
        radii_triso,
        radii_pebble,
        path,
        fuel_material,
        buffer_material,
        iPyC_material,
        SiC_material,
        oPyC_material,
        matrix_material,
        triso_PF=triso_PF,
        dist_triso_ini = dist_triso_ini
    )
    with open(path + "/geometry", "w") as f:
        f.write(block_plot)
        f.write(block_geom)
        f.write(block_pbed)
        f.write(block_triso)
    return True


if __name__ == "__main__":
    # Input
    Rin = 50
    Rout = 120
    H = 100
    Zmin = 0
    refl_thickness = 20
    shape = "cyl"
    name_pos = "fpb_pos"
    radii_pebble = [1.38, 1.80, 2.00]
    universe_pebbles = "u_pebbles"
    radii_triso = [0.02125, 0.03125, 0.03525, 0.03875, 0.04275]
    triso_PF = 0.22

    written_geometry = write_geometry_file(
        "./tmp",
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
    )
