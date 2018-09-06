

def perform(dirs, files, par, organizer):

    veBX, veBY = organizer.load_experiment_points()

    bX, bY = organizer.load_mesh_boundary_coordinates()

    baseData = organizer.load_red_file(files.virtual_experiment_data_files[0])

    from scipy.spatial import cKDTree
    import numpy as np

    fineMesh = cKDTree(baseData.data[:,:2])

    bIds = fineMesh.query(np.array([veBX, veBY]).T)[1]

    fineBoundaryMesh = cKDTree(np.array([veBX, veBY]).T)

    coarseInFine = fineBoundaryMesh.query(np.array([bX, bY]).T)[1]

    cBIds = bIds[coarseInFine]

    organizer.save(files.virtual_experiment_probe_mesh_ids, cBIds)

