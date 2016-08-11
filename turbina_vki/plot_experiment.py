import numpy as np
import os
import matplotlib.pyplot as plt
from plot_computation import provideVariableData


def expDataProvider(dirs, files, exp_variable):
    eXY = np.load(dirs.experiment+os.sep+"points.npy")
    for fname in files.experiment_datas:
        exp_data = np.ravel(np.load(dirs.experiment+os.sep+fname)[:, exp_variable])
        yield eXY[:, 0], exp_data


def scaleTo1(data):
    m = min(data)
    return (data - m) / (max(data) - m)

def perform(dirs, files, par):
    plt.figure()
    # for x, y in expDataProvider(dirs, files, 2):
    #     if max(y) < 1:
    #         y = scaleTo1(y)
    #         plt.plot(x, y, ".")

    for data in provideVariableData(files, par, varSelector=lambda d: d.ndata.mach, maxNum=1, minMach=0, maxMach=0.6, minAng=90, maxAng=94):
        # x = (data[0] - min(data[0]))*(max(data[0]) - min(data[0]))
        # y = scaleTo1(data[2])
        x = data[1]
        y = data[2]
        plt.plot(x, y, ".b")
        plt.plot(data[0],data[1],".r")
        print min(y)


    plt.show()

if __name__ == "__main__":
    import settings
    perform(settings.dirs, settings.files, settings.par)