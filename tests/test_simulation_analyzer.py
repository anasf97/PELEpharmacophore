import os
import shutil
import pytest

DIR = os.path.dirname(__file__)
SIMULATION_1 = os.path.join(DIR, "data/simulation_1")
TRAJECTORY_1 = os.path.join(DIR, "data/simulation_1/output/0/trajectory_1.pdb")


EXPECTED_IDS = ['NC1', 'NB1', 'NC3', 'O2', 'FD3', 'C1', 'CA5', 'CD1']

def test_get_grid_atoms(simulation_analyzer_1, atom_ids, trajectory=TRAJECTORY_1, expected_ids=EXPECTED_IDS):
    sa = simulation_analyzer_1
    str = sa.get_structure(trajectory)

    grid_atoms = sa.get_grid_atoms(str[0])

    assert atom_ids(grid_atoms) == expected_ids



EXPECTED_VOXELS_CV = [77, 582, 932, 1070, 1295, 1868, 2400, 2639]

def test_check_voxels(simulation_analyzer_1, active_voxels, trajectory=TRAJECTORY_1, expected_voxels=EXPECTED_VOXELS_CV):
    sa = simulation_analyzer_1
    str = sa.get_structure(trajectory)

    grid_atoms = sa.get_grid_atoms(str[0])
    grid = sa.check_voxels(grid_atoms)

    assert active_voxels(grid) == expected_voxels



EXPECTED_VOXELS_MG = [77, 582, 720, 932, 1070, 1267, 1295, 1646, 1868, 2400, 2414, 2626, 2639]

def test_merge_grids(simulation_analyzer_1, active_voxels, trajectory=TRAJECTORY_1, expected_voxels=EXPECTED_VOXELS_MG):
    sa = simulation_analyzer_1
    str = sa.get_structure(trajectory)

    grid_atoms1 = sa.get_grid_atoms(str[0])
    grid_atoms2 = sa.get_grid_atoms(str[1])

    grid1 = sa.check_voxels(grid_atoms1)
    grid2 = sa.check_voxels(grid_atoms2)

    final_grid = sa.merge_grids(grid1, grid2)

    assert active_voxels(final_grid) == expected_voxels



EXPECTED_VOXELS_AT = [77, 582, 720, 932, 1070, 1267, 1295, 1646, 1868, 2400, 2414, 2626, 2639]

def test_analyze_trajectory(simulation_analyzer_1, active_voxels, expected_voxels=EXPECTED_VOXELS_AT):
    sa = simulation_analyzer_1

    traj_grid = sa.analyze_trajectory(sa.traj_and_reports[0])

    assert active_voxels(traj_grid) == expected_voxels



EXPECTED_DICT = {'HBA': 1.1, 'ALI': 1.1, 'ARO': 0.6, 'HBD': 0.6}

def test_set_frequency_filter(simulation_analyzer_1, active_voxels, expected_dict=EXPECTED_DICT):
    sa = simulation_analyzer_1
    sa.grid = sa.analyze_trajectory(sa.traj_and_reports[0])

    sa.set_frequency_filter(1)

    assert sa.threshold_dict == expected_dict


EXPECTED_PHARMACOPHORES_PATH = os.path.join(DIR, "data/pharmacophores_1")

def test_save_pharmacophores(simulation_analyzer_1, compare_files, expected_pharmacophores_path=EXPECTED_PHARMACOPHORES_PATH):
    sa = simulation_analyzer_1
    sa.run(1)
    sa.set_frequency_filter(0)

    outdir = os.path.join(DIR, "result")
    sa.save_pharmacophores(outdir)

    result = [os.path.join(outdir, f) for f in os.listdir(outdir)]
    expected = [os.path.join(expected_pharmacophores_path, f) for f in os.listdir(expected_pharmacophores_path)]

    [compare_files(r, e) for r, e in zip(result, expected)]
    shutil.rmtree(outdir)
