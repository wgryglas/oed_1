__author__ = 'wgryglas'


def perform(dirs, files, par, organizer):
    """
    Script for creating symbolic links to each result (different angle of attack and outlet mach number) to store them
    in one place

    :param dirs: object containing directories information, comes from settings.py file
    :param files: object containing directories information, comes from settings.py file
    :param par: object containing directories information, comes from settings.py file
    :return: None
    """
    import os
    from wg.tools.system import ask_to_clean_dir

    ask_to_clean_dir(dirs.all_data)

    for d in os.listdir(dirs.converged_data):
        name = os.path.basename(d).replace('.', '-')
        src = dirs.converged_data+os.sep+d+os.sep+par.red_output_file_name
        out = dirs.all_data+os.sep+name+".dat"
        print "making symlink:", src, "->", out
        os.symlink(src, out)


if __name__ == "__main__":
    import settings
    perform(settings.dirs, settings.files, settings.par)