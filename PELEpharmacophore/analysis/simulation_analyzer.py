import os
import abc
import re
import glob
import numpy as np
from collections import OrderedDict
from itertools import accumulate, chain
import PELEpharmacophore.helpers as hl
import PELEpharmacophore.data.fragment_features as ff


class SimulationAnalyzer(metaclass=abc.ABCMeta):
    """
    Class for analysing PELE simulations.
    """

    def __init__(self, indir, features=None):
        """
        Create a new SimulationAnalyzer object.

        Parameters
        ----------
        indir : str
             Name of the simulation directory.
        """
        subdirs = [os.path.join(indir, subdir) for subdir in os.listdir(indir)]
        output_dir = os.path.join(indir, "output")

        if output_dir in subdirs:
            self.simulations = [Simulation(indir, features)]
        else:
            self.simulations = [Simulation(s, features) for s in subdirs]


    def set_ligand(self, chain, resname, resnum):
        """
        Set the parameters that define the ligand.

        Parameters
        ----------
        chain : str
             Ligand chain name.
        resname : str
             Ligand residue name.
        resnum : int
             Ligand residue number.
        """
        self.chain = chain
        self.resname = resname
        self.resnum = resnum


    def get_topology(self, file):
        """
        Parses a PDB file and returns a structure object.

        Parameters
        ----------
        file : str
            PDB file path.

        Returns
        ----------
        structure : Bio.PDB.Structure
            Biopython structure object.
        """
        return hl.load_topology(file)



    def get_indices(self, topology, resname, atomlist, first_index):
        """
        Gets all atoms defined in the `features` attribute.

        Parameters
        ----------
        model : Bio.PDB.Model
            Biopython model object.

        Returns
        ----------
        featured_grid_atoms : list of Atom objects
        """
        indlst  = [hl.get_indices(topology, resname, a) for a in atomlist]
        indices = np.concatenate([list(i) for i in indlst])
        res_indices =[i-first_index for i in indices]
        lengths = np.array([len(i) for i in indlst])

        return (res_indices , lengths)


    def get_coords(self, ncpus, steps):

        import tracemalloc

        tracemalloc.start()

        all_coord_dicts = []
        final_coord_dict={}

        for simulation in self.simulations:
          
            topology = self.get_topology(simulation.topfile)
            
            res_indices = topology.select(f"resname {self.resname}")

            first_index = res_indices[0]

            indices_dict = OrderedDict([(feature, self.get_indices(topology, self.resname, atomlist, first_index)) for feature, atomlist in simulation.features.items()])

            coord_dicts = hl.parallelize(get_coordinates, simulation.traj_and_reports, ncpus, indices_dict=indices_dict, resname=self.resname)


            sim_coord_dict = hl.merge_array_dicts(*coord_dicts)

            first_size, first_peak = tracemalloc.get_traced_memory()
            print(f"Loop memory usage is {first_size / 10**6}MB; Peak was {first_peak / 10**6}MB")

            #all_coord_dicts.append(sim_coord_dict)

            final_coord_dict = hl.merge_array_dicts(final_coord_dict, sim_coord_dict)

            first_size, first_peak = tracemalloc.get_traced_memory()
            print(f"Loop2 memory usage is {first_size / 10**6}MB; Peak was {first_peak / 10**6}MB")

        first_size, first_peak = tracemalloc.get_traced_memory()

        print(f"First memory usage is {first_size / 10**6}MB; Peak was {first_peak / 10**6}MB")

        #final_coord_dict = hl.merge_array_dicts(*all_coord_dicts)

        second_size, second_peak = tracemalloc.get_traced_memory()

        print(f"Second memory usage is {second_size / 10**6}MB; Peak was {second_peak / 10**6}MB")
        return final_coord_dict


    @abc.abstractmethod
    def save_pharmacophores(self):
        pass


def get_coordinates(traj_and_report, indices_dict, resname, steps):
    trajfile, report = traj_and_report
    indices = np.concatenate([i[0] for i in indices_dict.values()])
    print(indices)
    temp = np.argsort(indices)
    order_indices = np.empty_like(temp)
    order_indices[temp] = np.arange(len(indices))
    print(order_indices)
    accepted_steps = hl.accepted_pele_steps(report)

    coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=indices)
    coords = coords[:, order_indices]
    coords = coords[accepted_steps] # duplicate rows when a step is rejected
    coords = coords[:steps]

    coord_dict = {}
    start = 0
    for feature, (indices, lengths) in indices_dict.items():
        stop = start + len(indices)
        feature_coords = coords[:, start:stop, :]
        start = stop
        feature_coords = calc_cycle_centroids(feature_coords, lengths)
        coord_dict[feature] = feature_coords.reshape(-1, 3)

    return coord_dict

def calc_cycle_centroids(coords, lengths):
    ind = np.where(lengths > 1)[0]
    if ind.size == 0:
        return coords

    acc = list(accumulate(lengths))
    for i in ind:
        start = 0 if i == 0 else acc[i-1]
        stop = acc[i]
        cycle_coords = coords[:, start:stop, :]
        centroid = hl.centroid(cycle_coords)
        coords[:, start, :] = centroid
        coords[:, start+1:stop, :] = np.NaN
    coords = coords[~np.all(np.isnan(coords), axis=2)]
    return coords


class Simulation():
    """docstring for Simulation."""

    def __init__(self, indir, features=None):
        if features is None:
            frag_regex = r".*(?P<frag>frag\d+$)"
            frag = re.match(frag_regex, indir)['frag']
            self.features = ff.fragment_features[frag]

        else:
            self.features = features

        self.output = f"{indir}/output/"
        self.topfile = os.path.join(self.output, "topologies", "topology_0.pdb")
        self.trajectories = glob.glob(os.path.join(self.output, "0",  "trajectory_*.pdb"))
        self.reports = glob.glob(os.path.join(self.output, "0", "report_*"))
        self.traj_and_reports = self.match_traj_and_report()

    def match_traj_and_report(self):
        """
        Match each trajectory with its respective report.
        """
        self.trajectories.sort()
        self.reports.sort()
        traj_and_reports = list(zip(self.trajectories, self.reports))
        return traj_and_reports

    def set_features(d, fragment_features=ff.fragment_features):
        print(type(d))
        frag_regex = r".*(?P<frag>frag\d+$)"
        frag = re.match(frag_regex, d)['frag']
        features = fragment_features[frag]
        return features
