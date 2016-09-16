def perform(dirs, files, par, organizer):
    """
    Function makes links only to those results which are being assumed to be converged. All links are then stored
    in separate directory "converged" next to all numerical results.

    :param dirs: object containing directories information, comes from settings.py file
    :param files: object containing directories information, comes from settings.py file
    :param par: object containing directories information, comes from settings.py file
    :return: None
    """
    import os
    import re
    from glob import glob
    from wg.tools import file
    from wg.tools import system

    system.ask_to_clean_dir(dirs.converged_data)

    reg = re.compile('(ResL2*\s=*\s)(\d+\.\d+e(\+|\-)\d+)')

    for dir in os.listdir(dirs.root):
        if "aoa" not in dir:
            continue

        fname = dirs.root+os.sep+dir+os.sep+"*.log"
        fnames = glob(fname)
        if len(fnames) > 0:
            for line in file.ropen(fnames[0]):
                if "ResL2" in line:
                    res = reg.search(line)
                    if res and len(res.groups()) == 3 and res.group(3) == '-':
                        src = dirs.root+os.sep+dir
                        out = dirs.converged_data+os.sep+os.path.basename(dir)
                        print "making sym link ", src, "->", out
                        os.symlink(src, out)
                        break
                    elif res and len(res.groups())==3:
                        print res.groups(), dir
                        break




if __name__ == "__main__":
    import settings
    perform(settings.dirs, settings.files, settings.par)