import os
import shutil
import pytest
import numpy as np
from numpy.testing import assert_array_equal
import PELEpharmacophore.analysis.grid_analyzer as ga
import PELEpharmacophore.analysis.meanshift_analyzer as ma
import PELEpharmacophore.analysis.simulation_analyzer as sa

DIR = os.path.dirname(__file__)
SIMULATION_1 = os.path.join(DIR, "data/simulation_1")
TRAJECTORY_1 = os.path.join(DIR, "data/simulation_1/output/0/trajectory_1.pdb")
REPORT_1 = os.path.join(DIR, "data/simulation_1/output/0/report_1")
TOPOLOGY = os.path.join(DIR, "data/simulation_1/output/topologies/topology_0.pdb")
FEATURES =  {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': [('CA1', 'CA4'), ('CD1', 'CD4'), ('CC4', 'CC5', 'CC2')]}

ATOMLIST = ['NB1', 'NC3', 'O2']
EXPECTED_INDICES = [9, 17, 0]
EXPECTED_LENGTHS = [1, 1, 1]

def test_get_indices(create_analyzer, top_file=TOPOLOGY, atomlist = ATOMLIST, expected_indices=EXPECTED_INDICES, expected_lengths=EXPECTED_LENGTHS):
    print("first test")
    grid_a = create_analyzer(ga.GridAnalyzer)
    topology = grid_a.get_topology(top_file)
    res_indices = topology.select(f"resname {grid_a.resname}")
    first_index = res_indices[0]

    indices, lengths = grid_a.get_indices(topology, grid_a.resname, atomlist, first_index)
    assert np.all(indices == expected_indices)
    assert np.all(lengths == expected_lengths)


EXPECTED_COORD_DICT = dict([
                        ('HBD', np.array([[ 1.3959999, 14.473    , 28.174   ],
                                         [ 1.3959999, 14.473    , 28.174    ],
                                         [ 1.244    , 14.578    , 29.013    ]], dtype=np.float32)),
                        ('HBA', np.array([[ 0.49899998,  9.416     , 28.724998  ],
                                         [ 3.414     , 15.216     , 28.589     ],
                                         [-0.923     , 21.511     , 27.99      ],
                                         [ 0.49899998,  9.416     , 28.724998  ],
                                         [ 3.414     , 15.216     , 28.589     ],
                                         [-0.923     , 21.511     , 27.99      ],
                                         [ 0.41900003,  9.484     , 28.675999  ],
                                         [ 3.297     , 15.334001  , 29.105     ],
                                         [ 1.114     , 22.277     , 29.715     ]], dtype=np.float32)),
                        ('ALI', np.array([[ 8.342   , 11.331   , 29.886   ],
                                         [ 1.666   , 22.194   , 28.48    ],
                                         [ 8.342   , 11.331   , 29.886   ],
                                         [ 1.666   , 22.194   , 28.48    ],
                                         [ 8.273   , 11.280001, 29.374   ],
                                         [-0.438   , 21.522   , 27.604   ]], dtype=np.float32)),
                        ('ARO', np.array([[ 1.363    , 18.2915   , 27.817    ],
                                         [ 5.9279995, 12.535    , 29.3025   ],
                                         [ 2.592    , 14.249001 , 28.473333 ],
                                         [ 1.363    , 18.2915   , 27.817    ],
                                         [ 5.9279995, 12.535    , 29.3025   ],
                                         [ 2.592    , 14.249001 , 28.473333 ],
                                         [ 1.018    , 18.369    , 29.260498 ],
                                         [ 5.8475003, 12.571501 , 29.1915   ],
                                         [ 2.4653332, 14.354668 , 29.036001 ]], dtype=np.float32))
])

def test_get_coordinates(create_analyzer, compare_dicts, top_file=TOPOLOGY, features=FEATURES, traj=TRAJECTORY_1, report=REPORT_1, expected_coord_dict=EXPECTED_COORD_DICT):
    grid_a = create_analyzer(ga.GridAnalyzer)
    topology = grid_a.get_topology(top_file)
    traj_and_report = (traj, report)
    res_indices = topology.select(f"resname {grid_a.resname}")
    first_index = res_indices[0]

    indices_dict = {feature: grid_a.get_indices(topology, grid_a.resname, atomlist, first_index) \
                    for feature, atomlist in features.items()}
    coord_dict = sa.get_coordinates(traj_and_report, indices_dict=indices_dict, resname=grid_a.resname)
    compare_dicts(coord_dict,  expected_coord_dict)


COORDS = np.array([[ 1.3959999, 14.473    , 28.174    ],
                   [ 40.39599 , 14.473    , 28.174    ],
                   [ 1.244    , 14.578    , 29.013    ],
                   [ 1.3959999, 14.473    , 28.174    ],
                   [ 1.244    , 14.578    , 29.013    ]])

EXPECTED_GRID_ATOMS = np.array([[ 1.3959999, 14.473    , 28.174    ],
                                [ 1.244    , 14.578    , 29.013    ],
                                [ 1.3959999, 14.473    , 28.174    ],
                                [ 1.244    , 14.578    , 29.013    ]])

def test_get_grid_atoms(create_analyzer, coords=COORDS, expected_grid_atoms=EXPECTED_GRID_ATOMS):
    grid_a = create_analyzer(ga.GridAnalyzer)
    grid_atoms = grid_a.get_grid_atoms(coords)
    assert_array_equal(grid_atoms, expected_grid_atoms)

COORDS = np.array([[ 1.3959999, 14.473    , 28.174    ],
                   [ 1.244    , 14.578    , 29.013    ],
                   [ 1.3959999, 14.473    , 28.174    ],
                   [ 1.244    , 14.578    , 29.013    ]])

EXPECTED_VOXELS_INDS = [1070, 1267, 1070, 1267]

def test_check_voxels(create_analyzer, coords=COORDS, expected_voxel_indices=EXPECTED_VOXELS_INDS):
    grid_a = create_analyzer(ga.GridAnalyzer)
    voxel_centers = np.array([v.center for v in grid_a.grid.voxels])
    voxel_indices = ga.check_voxels(coords, voxel_centers)
    assert_array_equal(voxel_indices, expected_voxel_indices)

INPUT_INDICES = [1070, 1267, 1070, 1267, 1267]

FEATURE = "ARO"

EXPECTED_VOXEL_RESULTS = [(1070, {'ARO': 2}), (1267, {'ARO': 3})]

def test_fill_grid(create_analyzer, active_voxels, feature=FEATURE, input_indices=INPUT_INDICES, expected_voxel_results=EXPECTED_VOXEL_RESULTS):
    grid_a = create_analyzer(ga.GridAnalyzer)
    grid_a.fill_grid(feature, input_indices)
    voxel_results = active_voxels(grid_a.grid)
    assert voxel_results == expected_voxel_results


EXPECTED_DICT = {'HBA': 2.0, 'ALI': 2.0, 'ARO': 2.0, 'HBD': 2.0}

FILTER = 0

def test_set_frequency_filter_grid(run_analyzer, filter = FILTER, expected_dict=EXPECTED_DICT):
    grid_a = run_analyzer(ga.GridAnalyzer)

    grid_a.set_frequency_filter(FILTER)

    assert grid_a.threshold_dict == expected_dict

EXPECTED_PHARMACOPHORES_PATH = os.path.join(DIR, "data/pharmacophores_1")

def test_save_pharmacophores_grid(run_analyzer, compare_files, expected_pharmacophores_path=EXPECTED_PHARMACOPHORES_PATH):
    sa = run_analyzer(ga.GridAnalyzer)
    sa.set_frequency_filter(0)

    outdir = os.path.join(DIR, "result")
    sa.write_pdb_pharmacophores(outdir)

    result = [os.path.join(outdir, f) for f in os.listdir(outdir)]
    expected = [os.path.join(expected_pharmacophores_path, f) for f in os.listdir(expected_pharmacophores_path)]

    [compare_files(r, e) for r, e in zip(result, expected)]
    shutil.rmtree(outdir)

EXPECTED_CLUSTER_CENTERS =  dict([
                                ('HBD', np.array([[ 1.335, 14.514, 28.509]], dtype=np.float32)),
                                ('HBA', np.array([[ 0.467,  9.443, 28.705],
                                                  [ 3.367, 15.263, 28.795],
                                                  [-0.923, 21.511, 27.99 ],
                                                  [ 1.114, 22.277, 29.715]], dtype=np.float32)),
                                ('ALI', np.array([[ 8.314, 11.310, 29.681],
                                                  [ 1.666, 22.194, 28.480],
                                                  [-0.438, 21.522, 27.604]], dtype=np.float32)),
                                ('ARO', np.array([[ 1.363, 18.291, 27.817],
                                                  [ 5.895, 12.549, 29.258],
                                                  [ 2.541, 14.291, 28.698],
                                                  [ 1.018, 18.369, 29.260]], dtype=np.float32))
                                ])

EXPECTED_CLUSTER_FREQS = dict([
                               ('HBD', [5]),
                               ('HBA', [5, 5, 3, 2]),
                               ('ALI', [5, 3, 2]),
                               ('ARO', [3, 5, 5, 2]),

                                ])

def test_analyze_trajectory_meanshift(run_analyzer, compare_dicts, expected_cluster_centers=EXPECTED_CLUSTER_CENTERS,  expected_cluster_freqs=EXPECTED_CLUSTER_FREQS):
    ms_a = run_analyzer(ma.MeanshiftAnalyzer)
    ms_a.set_frequency_filter(0)

    cluster_centers, cluster_freqs = dict(), dict()
    for feature, clusters in ms_a.cluster_dict.items():
        centers = [cluster.center for cluster in clusters]
        freqs = [cluster.frequency for cluster in clusters]
        cluster_centers[feature] = centers
        cluster_freqs[feature] = freqs

    compare_dicts(cluster_centers, expected_cluster_centers)
    compare_dicts(cluster_freqs, expected_cluster_freqs)
