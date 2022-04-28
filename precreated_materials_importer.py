import shutil


def import_precreated(precreated_materials_path, dst, check=[]):
    shutil.copyfile(
        precreated_materials_path, dst+'/materials', follow_symlinks=True
    )
    if len(check) > 0:
        with open(precreated_materials_path) as f:
            lines = f.readlines()
        not_found = list(check)
        for line in lines:
            if line[:4] == "mat ":
                name_mat = line.split(" ")[1]
                if name_mat in check:
                    not_found.remove(name_mat)
        if len(not_found) > 0:
            raise Exception(
                "Precreated materials do not contain all necessary materials. Missing: {}".format(
                    not_found
                )
            )
    return True
