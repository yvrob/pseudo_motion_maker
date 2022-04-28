#%%
import math
import numpy as np
import struct
import os
import concurrent.futures
import shutil
from burnup_calculator import BU_to_compo_interpolator, calc_initial_BU
from glob import glob
import pandas as pd
from plotter import plot2D
from time import time
from itertools import repeat
from process_restart import *
import pickle


def PolarAngle(position):
    x = position[0]
    y = position[1]
    if x == 0.0 and y > 0.0:
        t = np.pi / 2.0
    elif x == 0.0 and y < 0.0:
        t = 3.0 * np.pi / 2.0
    elif x > 0.0 and y == 0.0:
        t = 0.0
    elif x < 0.0 and y == 0.0:
        t = np.pi
    elif x > 0.0 and y > 0.0:
        t = math.atan(y / x)
    elif x < 0.0 and y > 0.0:
        t = math.atan(y / x) + np.pi
    elif x < 0.0 and y < 0.0:
        t = math.atan(y / x) + np.pi
    elif x > 0.0 and y < 0.0:
        t = math.atan(y / x) + 2.0 * np.pi
    else:
        return 0.0
    return t


def replicate_dd(df, mpitasks):
    print("Replicating Serpent domain decomposition")
    positions = np.array(df[["x", "y", "z"]])

    subn = 100000
    nr = np.zeros((subn), int)
    nmat = len(positions)

    # Assign angular segment
    for i in positions:
        f = PolarAngle(i) / (2 * np.pi)
        f_truncated = f - float(int(f))
        n = int((float(subn) - 1.0) * f_truncated)
        nr[n] += 1

    # Assign MPI ID to segments
    MPI_ID = 0
    s = 0
    for n in range(subn):
        s += nr[n]
        nr[n] = MPI_ID
        if s > int(float(nmat) / float(mpitasks)) - 1:
            s = 0
            if MPI_ID < mpitasks - 1:
                MPI_ID += 1

    # Assign MPI ID to materials
    MPI_ID_list = []
    domains_indices = [[] for i in range(mpitasks)]
    for i, pos in enumerate(positions):
        f = PolarAngle(pos) / (2 * np.pi)
        f_truncated = f - float(int(f))
        n = int((float(subn) - 1.0) * f_truncated)
        if nr[n] < 0 or nr[n] > mpitasks - 1:
            raise Exception("indexing error: {} {}".format(nr[n], n))
        else:
            MPI_ID_list.append(nr[n])
            domains_indices[nr[n]].append(i)

    # Fill domains and check number of materials per domain
    s = 0
    MPI_cnt = [0 for i in range(mpitasks)]
    for MPI_ID in range(mpitasks):
        for i in MPI_ID_list:
            if i == MPI_ID:
                MPI_cnt[MPI_ID] += 1
        if MPI_cnt[MPI_ID] == 0:
            raise Exception(
                "Decomposition failed: no materials assigned to task {}".format(
                    MPI_ID + 1
                )
            )
        s += MPI_cnt[MPI_ID]
        print(
            "\tDomain {}: {} materials ({:.3f}%)".format(
                MPI_ID + 1, MPI_cnt[MPI_ID], 100.0 * (MPI_cnt[MPI_ID] / nmat)
            )
        )

    if s < nmat:
        print(0, "{} materials not decomposed".format(nmat - s))
    elif s > nmat:
        raise Exception("Not possible")

    # Assign to pebble
    df["MPI_ID"] = MPI_ID_list
    print("\tDone")
    return df


def write_domain_restart(df, MPI_ID, fuel_material, dep_output, files_out):
    pc = 0
    t0 = time()
    # print('Domain {}'.format(MPI_ID))
    df = df[df["MPI_ID"] == MPI_ID].sort_index(ascending=False)
    if dep_output[-2:] == ".m":
        interpolator, list_iso = BU_to_compo_interpolator(dep_output)
    else:
        with open(dep_output, "rb") as f:
            interpolator, list_iso = pickle.load(f)
    with open(files_out[MPI_ID], "wb") as fo:
        # Write actual materials
        cnt = 0
        pc_inc = 2

        # Parent
        name = "original"
        adens = interpolator(0)
        nnuc = adens.shape[0] - 1
        content = b""
        content += struct.pack("q", len(name))
        content += struct.pack("{}s".format(len(name)), str.encode(name))
        content += struct.pack("d", 0.0)  # BU global
        content += struct.pack("d", 0.0)  # BU days
        content += struct.pack("q", nnuc)
        content += struct.pack("d", adens[-1])
        content += struct.pack("d", 0.0)  # mdens
        content += struct.pack("d", 0.0)  # BU
        content += struct.pack("q", 0)  # pass
        content += struct.pack("q", -1)
        content += struct.pack("d", adens[-2])

        for j in range(len(adens) - 2):
            content += struct.pack("q", list_iso[j])
            content += struct.pack("d", adens[j])
        fo.write(content)

        # Zones
        for i in range(len(df)):
            mat = df.iloc[i]
            name = "{}z{}".format(fuel_material, mat.name + 1)
            adens = interpolator(mat["BU"])
            nnuc = adens.shape[0] - 1
            content = b""
            content += struct.pack("q", len(name))
            content += struct.pack("{}s".format(len(name)), str.encode(name))
            content += struct.pack("d", 0.0)  # BU global
            content += struct.pack("d", 0.0)  # BU days
            content += struct.pack("q", nnuc)
            content += struct.pack("d", adens[-1])
            content += struct.pack("d", 0.0)  # mdens
            content += struct.pack("d", mat["BU"])
            content += struct.pack("q", mat["passes"].astype(int))
            content += struct.pack("q", -1)
            content += struct.pack("d", adens[-2])

            for j in range(len(adens) - 2):
                content += struct.pack("q", list_iso[j])
                content += struct.pack("d", adens[j])
            fo.write(content)

            # Timer
            current_pc = (cnt + 1) / len(df) * 100
            if current_pc >= pc:
                if pc == 0:
                    print("\t[{}]\t{}%".format(MPI_ID, pc))
                else:
                    dt = int(time() - t0)
                    estimation = int(dt / (cnt + 1) * len(df))
                    print(
                        "\t[{}]\t{}%. Elapsed time: {}s. Estimation: {}s. Remaining time: {}s".format(
                            MPI_ID, pc, dt, estimation, estimation - dt
                        )
                    )
                while pc <= current_pc:
                    pc += pc_inc

            cnt += 1
        dt = int(time() - t0)
        print("\t[{}]\t100%. Elapsed time: {}s.".format(MPI_ID, dt))


def create_decomposed_restart_files(
    df,
    mpitasks,
    parallel,
    dep_output,
    path="./",
    name_out="first_compos",
    fuel_material="fuel",
):
    path = path + "/restart"
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

    files_out = []
    for i in range(mpitasks):
        files_out.append("{}/{}.wrk_dd{}".format(path, name_out, i))
        try:
            os.remove(files_out[-1])
        except:
            pass

    print("Creating decomposed restart files")
    if parallel:
        MPItasks = [i for i in range(mpitasks)]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for MPI_ID in zip(
                MPItasks,
                executor.map(
                    write_domain_restart,
                    repeat(df),
                    MPItasks,
                    repeat(fuel_material),
                    repeat(dep_output),
                    repeat(files_out),
                ),
            ):
                pass
    else:
        for MPI_ID in range(mpitasks):
            write_domain_restart(df, MPI_ID, fuel_material, dep_output, files_out)


if __name__ == "__main__":
    parallel = True
    mpitasks = 8
    fuel_material = "fuel_1"
    dep_output = "input_one_pebble_dep.m"
    name_out = "first_compos"

    #%% Clean/create folder
    average_discharge_burnup = 92  # MWd/kgHM
    npasses = 10
    residence_time = 300  # days
    max_bu_step = 10  # days
    H = 1100
    path = "./tmp/model"

    files = glob(path + "/" + "fpb_pos")
    files.sort(key=os.path.getmtime)
    df = pd.read_csv(
        files[-1], sep="\t", header=None, names=["x", "y", "z", "r_pebbles", "uni"]
    )
    df["dist"] = np.linalg.norm(df[["x", "y"]], axis=1)

    df = calc_initial_BU(df, residence_time, npasses, average_discharge_burnup, H)
    # Nsample = 10000
    # df = df.sample(n=Nsample, random_state=1)

    #%%%%%%%%%%%%

    df = replicate_dd(df, mpitasks)
    plot2D(df, dir_id=2, val=np.min(np.abs(df["z"])), field="MPI_ID")

    create_decomposed_restart_files(df, mpitasks, parallel, dep_output, path, name_out)
    #%%
    for i in range(mpitasks):
        file_in1 = glob(path + "/restart/*.wrk*")[i]
        restart1 = Restart_File(path_to_file=file_in1)
        restart1.read_restart(passes=True)
        restart1.write_text()
#%% Create new compositions
