"""
Creates fcc, hcp or bcc lattice filling dist_pebbles cuboid or (hollow) cylindrical volume based on dist_pebbles target packing fraction
"""

#%% Module
import numpy as np
import pandas as pd
import warnings
from plotter import *

#%% Functions
# Create lattice with correct lattice
def unit_cell(lattice_type, r_pebbles, nx, ny, nz):
    dist_pebbles = 2 * np.sqrt(2) * r_pebbles
    k, j, i = [
        v.flatten()
        for v in np.meshgrid(*([range(nz), range(ny), range(nx)]), indexing="ij")
    ]
    if lattice_type == "fcc":
        df = pd.DataFrame(
            {
                "x": dist_pebbles * i
                + dist_pebbles / 2 * (j % 2)
                + dist_pebbles / 2 * (k % 2),
                "y": dist_pebbles / 2 * j,
                "z": dist_pebbles / 2 * k,
            }
        )
    # elif lattice_type == 'hcp': # problem with nx,ny,nz to set better
    #     df = pd.DataFrame({
    #         'x': 2 * i + (j + k) % 2,
    #         'y': np.sqrt(3) * (j + 1/3 * (k % 2)),
    #         'z': 2 * np.sqrt(6) / 3 * k,
    #     })
    elif lattice_type == "sc":
        df = pd.DataFrame(
            {
                "x": 2 * r_pebbles * i,
                "y": 2 * r_pebbles * j,
                "z": 2 * r_pebbles * k,
            }
        )
    return df


# Fill volume with lattice and shpae
def fill_cube(
    lattice_type, dist_pebbles, sizeX, sizeY, sizeZ, centeredR=True, centeredZ=False
):
    nx = int(np.ceil(sizeX / (np.sqrt(2) * dist_pebbles)))
    ny = int(np.ceil(sizeY / (np.sqrt(2) * dist_pebbles)))
    nz = int(np.ceil(sizeZ / (np.sqrt(2) * dist_pebbles)))

    df = unit_cell(lattice_type, dist_pebbles, nx, ny, nz)
    if centeredR:
        df.x -= (df.x.max() - df.x.min()) / 2
        df.y -= (df.y.max() - df.y.min()) / 2
    if centeredZ:
        df.z -= (df.z.max() - df.z.min()) / 2
    return df


def shape_lattice(
    df, shape, r_pebbles, Rout, H, Rin=0, Zmin=0, *, partial_pebbles=False
):
    df["Rdist"] = np.linalg.norm(df[["x", "y"]], axis=1)
    df["dist"] = np.linalg.norm(
        df[["x", "y", "z"]] - df[["x", "y", "z"]].mean(), axis=1
    )
    if partial_pebbles:
        r_pebbles = 0
    if shape == "cyl":
        field = "Rdist"
    elif shape == "sph":
        field = "dist"
    if Rin > 0:
        df = df[df[field] >= Rin + r_pebbles]
    df = df[df[field] <= Rout - r_pebbles]
    if shape == "cyl":
        df = df[np.logical_and(df.z >= r_pebbles, df.z <= Zmin + H - r_pebbles)]
    df = df.reset_index().drop(["index"], axis=1)
    return df


def create_lattice(
    lattice_type,
    r_pebbles,
    dist_pebbles,
    shape,
    Rout,
    H,
    Rin=0,
    Zmin=0,
    *,
    partial_pebbles=False,
    centeredR=True,
    centeredZ=False,
    same_rows=True
):
    if shape == "cyl":
        sizeX = 2 * Rout
        sizeY = 2 * Rout
        sizeZ = H
    elif shape == "sph":
        sizeX = 2 * Rout
        sizeY = 2 * Rout
        sizeZ = 2 * Rout
    else:
        sizeX = Rout - Rin
        sizeY = Rout - Rin
        sizeZ = H

    # Fill cuboid
    df = fill_cube(
        lattice_type,
        dist_pebbles,
        sizeX,
        sizeY,
        sizeZ,
        centeredR=centeredR,
        centeredZ=centeredZ,
    )
    df["r_pebbles"] = r_pebbles

    # Shape to cylinder
    if shape == "cyl" or shape == "sph":
        df = shape_lattice(
            df, shape, r_pebbles, Rout, H, Rin, Zmin, partial_pebbles=partial_pebbles
        )
    if not isinstance(Zmin, type(None)) and Zmin != 0:
        df.z += Zmin
    shaped = False
    while not shaped:
        # Ensure same number of rows in fcc
        if same_rows and shape in ["cyl", "sq"] and lattice_type == "fcc":
            unique_z = np.unique(df.z)
            if len(unique_z) % 2 != 0:
                df = df[df.z < unique_z[-1]]

        if shape in ["cyl", "sq"]:
            lo_dZ = df.z.min() - Zmin
            hi_dZ = Zmin + H - df.z.max()
            dZ = (hi_dZ - lo_dZ) / 2
            df.z += dZ

        # TO DO BETTER
        if shape == "cyl":
            df = df[np.logical_and(df.z >= r_pebbles, df.z <= Zmin + H - r_pebbles)]

        if shape in ["cyl", "sq"]:
            lo_dZ = df.z.min() - Zmin
            hi_dZ = Zmin + H - df.z.max()
            dZ = (hi_dZ - lo_dZ) / 2
            df.z += dZ

        if shape == "cyl":
            df = df[np.logical_and(df.z >= r_pebbles, df.z <= Zmin + H - r_pebbles)]

        if same_rows and shape in ["cyl", "sq"] and lattice_type == "fcc":
            unique_z = np.unique(df.z)
            if len(unique_z) % 2 == 0:
                shaped = True
        else:
            shaped = True

    PF = calculate_packing_fraction(df, shape, r_pebbles, Rin, Rout, H)
    return df, PF


def calculate_packing_fraction(df, shape, r_pebbles, Rin, Rout, H):
    Vpebbles = len(df) * calculate_volume("sph", 0, r_pebbles, 0)
    Vcore = calculate_volume(shape, Rin, Rout, H)
    PF = Vpebbles / Vcore
    return PF


def calculate_volume(shape, Rin, Rout, H):
    if shape == "cyl":
        return np.pi * (Rout**2 - Rin**2) * H
    elif shape == "sph":
        return 4 / 3 * np.pi * (Rout**3 - Rin**3)
    else:
        return (Rout - Rin) ** 2 * H


# Iteratively finds right distance to match packing fraction
def find_lattice(
    lattice_type,
    r_pebbles,
    shape,
    Rout,
    H,
    target_PF,
    Rin=0,
    Zmin=0,
    dist_pebbles_ini=1000,
    mult_a=0.99,
    eps=1e-2,
    *,
    partial_pebbles=False,
    centeredR=True,
    centeredZ=False,
    same_rows=True,
    verbose=False
):
    # Quick check
    if (
        lattice_type == "fcc"
        and target_PF > np.pi * np.sqrt(2) / 6
        or lattice_type == "bc"
        and target_PF > np.pi / 6
        or lattice_type == "hcp"
        and target_PF > np.pi * np.sqrt(2) / 6
    ):
        raise Exception(
            "Cannot obtain such dist_pebbles high packing fraction with this lattice type: {}".format(
                lattice_type
            )
        )

    # Initializers
    PF = 0
    found = False
    trying_again = False
    dist_pebbles = dist_pebbles_ini
    cnt = 0
    # Iterate (TO CHECK)
    while not found and mult_a < 1:
        while PF < target_PF or trying_again:
            trying_again = False
            df, PF = create_lattice(
                lattice_type,
                r_pebbles,
                dist_pebbles,
                shape,
                Rout,
                H,
                Rin,
                Zmin,
                partial_pebbles=partial_pebbles,
                centeredR=centeredR,
                centeredZ=centeredZ,
                same_rows=same_rows,
            )
            if abs(PF - target_PF) / target_PF < eps:
                found = True
                print(
                    "Found! dist_pebbles={:.10f}, PF={:.4f}, N={}".format(
                        dist_pebbles, PF, len(df)
                    )
                )
                break
            dist_pebbles *= mult_a
        if not found:
            if verbose:
                print(
                    "Not found with dist_pebbles={:.10f}, mult_a={:.4f} (PF={:.4f}, N={}), trying again from dist_pebbles={:.4f} with mult_a={:.4f}".format(
                        dist_pebbles / mult_a,
                        mult_a,
                        PF,
                        len(df),
                        dist_pebbles / mult_a**2,
                        mult_a + (1 - mult_a) / 2,
                    )
                )
            dist_pebbles /= mult_a**2
            mult_a += (1 - mult_a) / 2
            trying_again = True
            cnt += 1
        if cnt > 50:
            break
    return df, dist_pebbles, PF, found


if __name__ == "__main__":
    #%% Script
    # Input parameters
    # r_pebbles = 2  # pebbles radius
    # Rin = 40  # inner radius, 0 if no inner radius
    # Rout = 120  # 120 # core radius
    # Zmin = 10  # Minimum elevation
    # H = 310  # 310 # core height
    # target_PF = 0.6  # packing fraction to reach
    # lattice_type = "fcc"  # can be hcp (not usable), fcc, sc
    # shape = "cyl"  # can be cyl or sph whatever (will be cuboid otherwise)

    r_pebbles = 4.275000e-02  # pebbles radius
    Rin = 1.38  # inner radius, 0 if no inner radius
    Rout = 1.80  # 120 # core radius
    Zmin = None  # Minimum elevation
    H = None  # 310 # core height
    target_PF = 0.22  # packing fraction to reach
    lattice_type = "sc"  # can be hcp (not usable), fcc, sc
    shape = "sph"  # can be cyl or sph whatever (will be cuboid otherwise)

    # Secondary input
    same_rows = True  # needed for discrete motion
    partial_pebbles = False  # allow for pebbles crossing bounds in volume
    centeredR = True  # To radially center the geometry
    centeredZ = False

    # Parameters for iterative search with target packing fraction within volume
    dist_pebbles_ini = r_pebbles * 5  # initial distance, huge on purpose
    mult_a = 0.9  # initial multiplier for dist_pebbles, everytime the PF is too small
    eps = 1e-3  # desired precision

    # Iterate
    df, dist_pebbles, PF, found = find_lattice(
        lattice_type,
        r_pebbles,
        shape,
        Rout,
        H,
        target_PF,
        Rin,
        Zmin,
        dist_pebbles_ini,
        mult_a,
        eps,
        partial_pebbles=partial_pebbles,
        centeredR=centeredR,
        centeredZ=centeredZ,
        same_rows=same_rows,
        verbose=True,
    )
    if not found:
        warnings.warn(
            "Did not find any good solution! Taking best with dist_pebbles={:.4f}, (PF={:.4f})".format(
                dist_pebbles, PF
            )
        )

    for i in range(3):
        direction = ["x", "y", "z"][i]
        value = np.array(df[direction]).flat[
            np.abs(df[direction] - df[direction].mean()).argmin()
        ]
        plot2D(df, i, value, "Rdist")

    plot_df(df, "Rdist")
