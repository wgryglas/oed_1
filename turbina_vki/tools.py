__author__ = 'wgryglas'


def __performe_base_experiment_transofrm__(X, Y):
    """
    Function which performs base transform of boundary coordinates
    measured in experiment to coordinates used in computation
    :param X: probe X coordinates
    :param Y: probe Y coordinates
    :return: x,y arrays of transformed coordinates
    """
    import numpy as np
    angle = np.radians(-26.7)
    c = np.cos(angle)
    s = np.sin(angle)
    dy = 0.534
    return X*c + Y*s, X*s - Y*c + dy


def fit_experiment_porbe_locations_to_mesh_boundary(probeX, probeY, boundaryX, boundaryY):
    """
    Function converting experiment probe locations into mesh boundary curve defined by boundary nodes
    :param probeX:
    :param probeY:
    :param boundaryX:
    :param boundaryY:
    :return: x,y coordinates of probe locations transformed onto boundary
    """
    from wg.tools.geom import fit_coords

    px, py = __performe_base_experiment_transofrm__(probeX, probeY)
    px, py, _ = fit_coords(px, py, boundaryX, boundaryY)
    return px, py
