__author__ = 'wgryglas'


def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED import POD
    import matplotlib.pyplot as plt
    eigens = organizer.load(files.eigen_values)

    plt.figure()
    POD.Modes.add_spectrum_plot(eigens, par.num_modes)
    plt.show()


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)