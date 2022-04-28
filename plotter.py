#%%
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits import mplot3d
from glob import glob
import os

#%%
# Plot
def slicing(df, dir_id, val):
    return df[np.abs(df[["x", "y", "z"][dir_id]] - val) == 0]


def plot2D(df, dir_id, val, field=None, equal=True):
    if dir_id == 0:
        xdir = 1
        ydir = 2
    elif dir_id == 1:
        xdir = 0
        ydir = 2
    elif dir_id == 2:
        xdir = 0
        ydir = 1

    df = slicing(df, dir_id, val)
    x = np.array(df[["x", "y", "z"][xdir]])
    y = np.array(df[["x", "y", "z"][ydir]])
    r_pebbles = np.array(df["r_pebbles"])

    patches = []
    for i in range(len(df)):
        circle = Circle((x[i], y[i]), r_pebbles[i])
        patches.append(circle)

    if isinstance(field, type(None)):
        colors = r_pebbles
    else:
        colors = np.array(df[field])

    p = PatchCollection(patches)
    p.set_array(colors)
    ax = plt.gca()
    ax.add_collection(p)
    plt.xlabel(["x", "y", "z"][xdir])
    plt.ylabel(["x", "y", "z"][ydir])
    plt.title("{}={:.3f}".format(["x", "y", "z"][dir_id], val))
    ax.autoscale_view()
    if equal:
        plt.gca().set_aspect("equal", adjustable="box")


def plot_df(df_to_plot, field, view=[-140, 60], scatter_size=10, alpha=1):
    fig = plt.figure(figsize=(10, 10))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()
    ax = fig.add_subplot(111, projection="3d")  # , proj_type = 'ortho')
    ax.scatter3D(
        df_to_plot.x,
        df_to_plot.y,
        df_to_plot.z,
        s=scatter_size,
        c=df_to_plot[field],
        alpha=alpha,
        zorder=1,
    )
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_zlabel("z [cm]")
    mult = 1.1
    ax.set_xlim3d(mult * min(df_to_plot.x), mult * max(df_to_plot.x))
    ax.set_ylim3d(mult * min(df_to_plot.y), mult * max(df_to_plot.y))
    ax.set_zlim3d(mult * min(df_to_plot.z), mult * max(df_to_plot.z))
    ax.view_init(view[0], view[1])


def import_last(path, pattern, plot=False):
    files = glob(path + "/" + pattern)
    files.sort(key=os.path.getmtime)
    df = pd.read_csv(
        files[-1], sep="\t", header=None, names=["x", "y", "z", "r_pebbles", "uni"]
    )
    df["dist"] = np.linalg.norm(df[["x", "y", "z"]], axis=1)
    if plot:
        plot_df(df, "dist", alpha=0.4)
        plt.title(pattern.replace("*", "").split(".")[0])
        plt.figure(figsize=(10, 3))
        for i in range(3):
            plt.subplot(1, 3, i + 1)
            direction = ["x", "y", "z"][i]
            value = np.array(df[direction]).flat[
                np.abs(df[direction] - df[direction].mean()).argmin()
            ]
            plot2D(df, i, value, "dist")
        plt.tight_layout()
        plt.show()
    return df
