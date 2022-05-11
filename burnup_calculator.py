import pandas as pd
import numpy as np
from glob import glob
from plotter import *
import pickle


def calculate_nrows_step(df, max_bu_step, residence_time, npasses):
    bu_step = np.inf
    nrows_step = len(np.unique(df.z))
    while bu_step > max_bu_step:
        bu_step = burnstep(df, residence_time, npasses, nrows_step)
        nrows_step -= 1
    return nrows_step


def burnstep(df, residence_time, npasses, nrows_step):
    df_sample = df[np.logical_and(df.x==df.x[0], df.y==df.y[0])]
    h = df_sample.z.max()-df_sample.z.min()
    dh = np.unique(df_sample.z)[1] - np.unique(df_sample.z)[0]
    t = residence_time/npasses
    dt = t*dh*nrows_step/(h+dh)
    return dt

def add_random_passes(df, npasses):
    passes = np.array_split([i for i in range(len(df))], npasses)
    for i in range(len(passes)):
        passes[i][:] = i + 1
    passes = np.concatenate(passes).ravel()
    np.random.shuffle(passes)
    df["passes"] = passes.astype(int)
    return df


def z_pass_to_time(df, residence_time, H, direction):
    npasses = df.passes.max()
    time_per_pass = residence_time / npasses
    previous_passes_time = time_per_pass * (df.passes - 1)
    if direction == +1:
        Zmin = df.z.min()
        current_pass_time = (df.z - Zmin) / H * time_per_pass
    elif direction == -1:
        Zmax = df.z.max()
        current_pass_time = (Zmax - df.z) / H * time_per_pass
    df["time"] = previous_passes_time + current_pass_time
    return df


def time_to_BU(df, residence_time, average_discharge_burnup):
    df["BU"] = df["time"] / residence_time * average_discharge_burnup
    return df


def calc_initial_BU(df, residence_time, npasses, average_discharge_burnup, H, direction):
    df = add_random_passes(df, npasses)
    df = z_pass_to_time(df, residence_time, H, direction)
    df = time_to_BU(df, residence_time, average_discharge_burnup)
    return df


def BU_to_compo_interpolator(dep_output):
    import serpentTools as st
    from scipy.interpolate import interp1d

    dep = st.read(dep_output)
    mat = dep.materials["fuel_1"]
    list_iso = mat.zai
    interpolator = interp1d(mat.burnup, mat.adens, kind="linear")
    with open("{}.interpolator".format(dep_output.split(".m")[0]), "wb") as file_handle:
        pickle.dump((interpolator, list_iso), file_handle)
    return interpolator, list_iso


if __name__ == "__main__":
    path = "./tmp/model"

    average_discharge_burnup = 92  # MWd/kgHM
    npasses = 10
    residence_time = 300  # days
    max_bu_step = 10  # days
    H = 1100

    files = glob(path + "/" + "fpb_pos")
    files.sort(key=os.path.getmtime)
    df = pd.read_csv(
        files[-1], sep="\t", header=None, names=["x", "y", "z", "r_pebbles", "uni"]
    )
    df["dist"] = np.linalg.norm(df[["x", "y"]], axis=1)

    df = calc_initial_BU(df, residence_time, npasses, average_discharge_burnup, H, direction)
    plot2D(df, 0, 5.68434e-14, field="BU", equal=False)
    df2 = (
        df[df.passes == 10]
        .groupby(
            pd.cut(
                df[df.passes == 10]["dist"],
                np.linspace(df.dist.min(), df.dist.max(), 30),
            )
        )["BU"]
        .mean()
    )
    plt.figure()
    df2.plot()
    df2 = (
        df[df.passes == 10]
        .groupby(
            pd.cut(df[df.passes == 10]["z"], np.linspace(df.z.min(), df.z.max(), 30))
        )["BU"]
        .mean()
    )
    plt.figure()
    df2.plot()
