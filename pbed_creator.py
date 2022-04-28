import numpy as np
import pandas as pd


def create_pbed(pos_list, r_list, uni_list):
    if isinstance(uni_list, str):
        universe = str(uni_list)
        uni_list = np.empty_like(r_list, dtype="U25")
        uni_list[:] = universe
    elif len(uni_list) == 1:
        universe = uni_list[0]
        uni_list = np.empty_like(r_list, dtype="U25")
        uni_list[:] = universe

    pos_list = pd.DataFrame(pos_list, columns=["x", "y", "z"])
    r_list = pd.DataFrame(r_list, columns=["r"])
    uni_list = pd.DataFrame(uni_list, columns=["uni"])
    df = pd.concat([pos_list, r_list, uni_list], axis=1)
    return df


def pbed_to_file(df, path, name, precision=5):
    df.to_csv(path + "/" + name, sep="\t", index=False, header=False, float_format='%.{}E'.format(precision))


def write_pbed(path, name, pos_list, r_list, uni_list):
    df = create_pbed(pos_list, r_list, uni_list)
    pbed_to_file(df, path, name)
    return True


if __name__ == "__main__":
    # Input
    universe = "fuel"
    path = "./tmp"
    name = "fpb_pos"

    pos = pd.DataFrame(
        np.hstack((np.random.random((100, 3)), np.ones((100, 1)) * 2)),
        columns=["x", "y", "z", "r"],
    )
    uni_list = np.empty_like(r_list, dtype="U25")
    uni_list[:] = universe

    write_pbed(path, name, pos[["x", "y", "z"]], pos["r"], uni_list)
