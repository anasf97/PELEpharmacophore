import os
import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import PELEpharmacophore.analysis.grid_analyzer as ga
import PELEpharmacophore.analysis.meanshift_analyzer as ma


DIR = os.path.dirname(__file__)
SIMULATION = os.path.join(DIR, "data/simulation_1")

CHAIN, RESNAME, RESNUM = "L", "SB2", 800

CENTER=[2.173, 15.561, 28.257]
RADIUS=7

FEATURES =  {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': [('CA1', 'CA4'), ('CD1', 'CD4'), ('CC4', 'CC5', 'CC2')]}

@pytest.fixture
def create_analyzer():
    def _create_analyzer(analyzer_class):
        a = analyzer_class(SIMULATION)
        a.set_ligand(CHAIN, RESNAME, RESNUM)
        a.set_features(FEATURES)
        if analyzer_class == ga.GridAnalyzer:
            a.set_grid(CENTER, RADIUS)
        return a
    return _create_analyzer

@pytest.fixture
def run_analyzer(create_analyzer):
    def _run_analyzer(analyzer):
        a = create_analyzer(analyzer)
        a.run(1)
        return a
    return _run_analyzer

@pytest.fixture
def active_voxels():
    def _active_voxels(grid):
        return [(i, v.freq_dict) for i, v in enumerate(grid.voxels) if v.freq_dict]
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

@pytest.fixture
def compare_dicts():
    def _compare_dicts(dict1, dict2):
        assert set(dict1.keys()) == set(dict2.keys())
        [assert_allclose(dict1[key], dict2[key], rtol=1e-3) for key in dict1]
    return _compare_dicts
