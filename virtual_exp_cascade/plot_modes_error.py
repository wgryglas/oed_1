
def perform(dirs, files, par, organizer):
    import numpy as np
    import itertools
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('text', usetex=True)

    modes = [m for m in organizer.load_modes()]
    #marker = itertools.cycle(('.', 'x', 'o', '+', ','))
    markers = list()

    # for m in Line2D.markers:
    #     try:
    #         if len(m) == 1 and m != ' ':
    #             markers.append(m)
    #     except TypeError:
    #         pass
    #
    # print markers

    legends = []
    lines = []

    fig = plt.figure(**organizer.get_figure_setup())

    for dfile in files.modes_test_data_file:
        legends.append(dfile)
        dataRed = organizer.load_red_file(files.modes_test_data_file[dfile])
        var = par.optimization_variable
        data = dataRed.ndata[var]
        dataRange = np.abs(max(data) - min(data))
        erros = list()


        for i in range(len(modes)):

            print "nmodes", i

            A = np.array([m.ndata[var] for m in modes[:i+1]]).T
            coeffs = np.linalg.lstsq(A, data)[0]
            recon = A.dot(coeffs)

            diff = max(np.abs(data - recon))/dataRange

            erros.append(diff)

        line,  = plt.plot(range(len(modes)), erros, "-"+organizer.next_marker, color="black", markersize=10, markerfacecolor='None')
        lines.append(line)

    plt.legend(lines, legends)
    plt.grid(True)

    plt.savefig(dirs.figures + "/modes_reconstruction_comparison.png")

