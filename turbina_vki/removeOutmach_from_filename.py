

def perform(dirs, files, par):

    from os import rename
    from os.path import join
    from glob import glob1

    for d in glob1(dirs.all_data, "*outmach*.dat"):
        nd = "-".join(d.split("-")[:-1])[:-8]+".dat"
        rename(join(dirs.all_data, d), join(dirs.all_data, nd))


# if __name__ == "__main__":
#     from settings import *
#     perform(dirs, files, par)