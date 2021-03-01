import os
import pytest
import PELEpharmacophore.analysis.simulation_analyzer as sa


DIR = os.path.dirname(__file__)
SIMULATION_1 = os.path.join(DIR, "data/simulation_1")
TRAJECTORY_1 = os.path.join(DIR, "data/simulation_1/output/0/trajectory_1.pdb")

CHAIN_1, RESNAME_1, RESNUM_1 = "L", "SB2", 800

CENTER_1=[2.173, 15.561, 28.257]
RADIUS_1=7

FEATURES_1 = {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': ['CA5', 'CD1']}


@pytest.fixture
def simulation_analyzer_1():
    s = sa.SimulationAnalyzer(SIMULATION_1)
    s.set_ligand(CHAIN_1, RESNAME_1, RESNUM_1)
    s.set_features(FEATURES_1)
    s.set_grid(CENTER_1, RADIUS_1)
    return s


@pytest.fixture
def atom_ids():
    def _atom_ids(atoms):
        return [a.atom.id for a in atoms]
    return _atom_ids


@pytest.fixture
def active_voxels():
    def _active_voxels(grid):
        return [i for i, v in enumerate(grid.voxels) if v.freq_dict]
    return _active_voxels


@pytest.fixture
def compare_files():
    def _compare_files(file1, file2):
        with open(file1, 'r') as f1:
            lines1 = [line for line in f1.readlines()]
        with open(file2, 'r') as f2:
            lines2 = [line for line in f2.readlines()]
        assert len(lines1) == len(lines2), \
            'Number of lines do not match: ' \
            + str(len(lines1)) + ' and ' + str(len(lines2))
        for i, (line1, line2) in enumerate(zip(lines1, lines2)):
            assert line1 == line2, \
                'Found different lines at line {}:'.format(i) + '\n' + line1 + line2
    return _compare_files
