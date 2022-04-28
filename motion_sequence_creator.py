#%%
from plotter import *
import shutil
from burnup_calculator import burnstep, calculate_nrows_step
from time import time
import matplotlib.cm as cm


def columns_rows(df):
    columns_group = df.groupby(["y", "x"], sort=False)
    columns_id = columns_group.ngroup().add(1) - 1
    df["column_id"] = columns_id

    rows_group = df.groupby("z")
    rows_id = rows_group.ngroup().add(1) - 1
    df["row_id"] = rows_id
    rows_group = df.groupby(df["row_id"] // 2)
    rows_id = rows_group.ngroup().add(1) - 1
    df["row_id"] = rows_id
    return df


def create_sequence_scheme(df, max_bu_step, residence_time, npasses, ncycles):
    nrows_step = calculate_nrows_step(df, max_bu_step, residence_time, npasses)
    df = columns_rows(df)
    nrows = df.row_id.max() + 1
    fraction = nrows_step / nrows
    npasses_tot = ncycles * npasses
    return [[npasses_tot, fraction, fraction]]


def calculate_steps_from_sequence(df, sequence_scheme):
    #%% Step calculation
    df = columns_rows(df)
    nrows = df.row_id.max() + 1
    steps = []
    for i in sequence_scheme:
        # print(i)
        s = 0
        while s < nrows * i[0]:
            step = int(nrows * np.random.uniform(i[1], i[2] + 1e-10))
            if step == 0:
                step = 1
            # steps.append(int(nrows*i[1]))
            steps.append(step)
            s += step
    return steps


def create_bustep(df, sequence_scheme, residence_time, npasses):
    steps = calculate_steps_from_sequence(df, sequence_scheme)
    calculated_burnups = [burnstep(df, residence_time, npasses, i) for i in steps]
    budays = np.concatenate(([0], np.cumsum(calculated_burnups)))
    s = "% Total time: {:.5} days (= {:.5} years)\n".format(
        sum(calculated_burnups), sum(calculated_burnups) / 365
    )
    s += "dep daystep {}\n".format(" ".join(np.array(calculated_burnups).astype(str)))
    return s


def write_bustep(path, df, sequence_scheme, residence_time, npasses, own_file=False):
    string = create_bustep(df, sequence_scheme, residence_time, npasses)
    if own_file:
        with open(path, "w") as f:
            f.write(string)
    else:
        with open(path, "a") as f:
            f.write(string + "\n")


def show_list_index(df, list_index, dir_id, value):
    fig, ax = plt.subplots()
    labels = np.array(["x", "y", "z"])
    direction = labels[dir_id]
    data = df[df[direction] == value]
    array = np.array(list_index)
    array = array[data.index]
    if direction == "x":
        indices = [1, 2]
    elif direction == "y":
        indices = [0, 2]
    elif direction == "z":
        indices = [0, 1]
    x = data[labels[indices[0]]]
    y = data[labels[indices[1]]]
    label = labels[indices]
    scat = ax.scatter(x, y, s=3, c=array, cmap=cm.jet)
    ax.set_xlabel("{} [cm]".format(label[0]))
    ax.set_ylabel("{} [cm]".format(label[1]))
    cbar = plt.colorbar(scat, shrink=0.8)
    plt.title("Material indices")
    ax.axis("equal")


def show_columns(df, list_index, direction, value):
    fig, ax = plt.subplots()
    labels = np.array(["x", "y", "z"])
    direction = labels[dir_id]
    data = df[df[direction] == value]
    array = np.array(df.reindex(np.array(list_index)).column_id)
    array = array[data.index]
    if direction == "x":
        indices = [1, 2]
    elif direction == "y":
        indices = [0, 2]
    elif direction == "z":
        indices = [0, 1]
    x = data[labels[indices[0]]]
    y = data[labels[indices[1]]]
    label = labels[indices]
    scat = ax.scatter(
        x, y, s=3, c=array, cmap=cm.jet, clim=[np.min(array), np.max(array)]
    )
    ax.set_xlabel("{} [cm]".format(label[0]))
    ax.set_ylabel("{} [cm]".format(label[1]))
    cbar = plt.colorbar(scat, shrink=0.8)
    plt.title("Column indices")
    ax.axis("equal")


def write_motion_sequence(
    path,
    df,
    direction,
    sequence_scheme,
    lattice_type,
    npasses,
    residence_time,
    randomize_reloading=True,
    plot=False,
):
    #%% Clean/create folder
    if os.path.exists(path + "/indices"):
        shutil.rmtree(path + "/indices")
    os.makedirs(path + "/indices")

    df = columns_rows(df)
    nrows = df.row_id.max() + 1
    ncols = df.column_id.max() + 1

    steps = calculate_steps_from_sequence(df, sequence_scheme)
    Nsteps = len(steps)

    #%% Create initial Matrix of pebbles based on columns/rows
    indices_matrix = np.ones((nrows, ncols), dtype=int) * -100000
    indices_matrix[df.row_id, df.column_id] = df.index

    indices_matrix_original = np.array(indices_matrix)
    list_index = indices_matrix[df.row_id, df.column_id]

    print("Starting moving pebbles")
    print("Step", 0)
    np.savetxt(
        path + "/indices/input_mat_indices_{}".format(0),
        list_index,
        fmt="%d",
        delimiter="\n",
    )
    with open(path + "/indices/input_mat_indices_recirculation_{}".format(0), "w"):
        pass
    if plot:
        show_list_index(df, list_index, 0, 0)
        show_columns(df, list_index, 0, 0)
        plt.show()
    t0 = time()
    for i_step in range(Nsteps):
        step = steps[i_step]
        if lattice_type == "sc":
            Nrows_to_move = step
        elif lattice_type == "fcc":
            Nrows_to_move = step

        dt = time() - t0
        estimation = dt / (i_step + 1) * Nsteps
        print(
            "Step {} / {}. Elapsed time: {:.1f}s. Estimated remaining time: {:.1f}s".format(
                i_step + 1, Nsteps, dt, estimation - dt
            )
        )
        # if i_step%10 == 0:
        #     clear_output(wait=True)

        #%% Motion
        indices_matrix = np.roll(indices_matrix, Nrows_to_move * direction, 0)
        ncolumns = indices_matrix.shape[1]

        #%% Determine recirculating rows
        if direction == +1:
            recirculated_rows = [i for i in range(0, min(Nrows_to_move, nrows))]
            recirculating_rows = [
                i for i in range(nrows - min(Nrows_to_move, nrows), nrows)
            ]
        elif direction == -1:
            recirculated_rows = [
                i for i in range(nrows - min(Nrows_to_move, nrows), nrows)
            ]
            recirculating_rows = [i for i in range(0, min(Nrows_to_move, nrows))]

        #%% Shuffle if necessary
        if randomize_reloading:
            print(
                "\tShuffling {} rows: {} ".format(
                    len(recirculated_rows), recirculated_rows
                )
            )
            for row in recirculated_rows:
                indices_matrix[row, :] = indices_matrix[
                    row, np.random.permutation(ncolumns)
                ]

        #%% Record new positions
        list_index = indices_matrix[df.row_id, df.column_id]

        if Nrows_to_move > 0:
            list_recirculation = np.concatenate(
                indices_matrix_original[recirculating_rows, :]
            ).ravel()
        else:
            list_recirculation = []

        np.savetxt(
            path + "/indices/input_mat_indices_{}".format(i_step + 1),
            list_index,
            fmt="%d",
            delimiter="\n",
        )
        np.savetxt(
            path + "/indices/input_mat_indices_recirculation_{}".format(i_step + 1),
            list_recirculation,
            fmt="%d",
            delimiter="\n",
        )
        if plot:
            show_list_index(df, list_index, 0, 0)
            show_columns(df, list_index, 0, 0)
            plt.show()


if __name__ == "__main__":
    path = "./tmp/model"
    df = import_last(path, pattern="fpb_pos")
    residence_time = 300
    npasses = 10
    lattice_type = "fcc"

    direction = -1  # +1=up, -1=down
    randomize_reloading = True
    sequence_scheme = [(24, 0.1, 0.1)]
    plot = True

    write_bustep("./test_bu", df, sequence_scheme, residence_time, npasses)
    write_motion_sequence(
        path,
        df,
        direction,
        sequence_scheme,
        lattice_type,
        npasses,
        residence_time,
        randomize_reloading=True,
        plot=plot,
    )

    # %%
    create_bustep(calculated_burnups)
