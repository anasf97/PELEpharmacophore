import os
import abc
import re
import glob
import numpy as np
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



    def get_indices(self, topology, resname, atomlist):
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
        lengths = np.array([len(i) for i in indlst])

        return (indices , lengths)


    @abc.abstractmethod
    def run(self, ncpus):
        import tracemalloc

        tracemalloc.start()
        all_coord_dicts = []

        for simulation in self.simulations:

            print(simulation.output)

            topology = self.get_topology(simulation.topfile)

            indices_dict = {feature: self.get_indices(topology, self.resname, atomlist) \
                            for feature, atomlist in simulation.features.items()}

            print(indices_dict)

            coord_dicts = hl.parallelize(get_coordinates, simulation.traj_and_reports, ncpus, indices_dict=indices_dict)
            
            print(coord_dicts)	    
            gen_dict = hl.gen_array_dicts(*coord_dicts)
	   
            print(gen_dict)

            #all_coord_dicts.append(gen_dict)

        first_size, first_peak = tracemalloc.get_traced_memory()

        print(f"Second memory usage is {first_size / 10**6}MB; Peak was {first_peak / 10**6}MB")

       # merged_all_coord_dict = hl.gen_array_dicts(*all_coord_dicts)

        #final_coord_dict = {feature: chain.from_iterable(gen) \
        #                    for feature, gen in merged_all_coord_dict.items()}

        #second_size, second_peak = tracemalloc.get_traced_memory()

        print(f"Second memory usage is {second_size / 10**6}MB; Peak was {second_peak / 10**6}MB")
        #return final_coord_dict
        return

    @abc.abstractmethod
    def save_pharmacophores(self):
        pass

def get_simulation_coordinates(simulation, resname):
        topology = self.get_topology(simulation.topfile)

        indices_dict = {feature: self.get_indices(topology, resname, atomlist) \
                        for feature, atomlist in simulation.features.items()}

        coord_dicts = hl.parallelize(get_coordinates, simulation.traj_and_reports, 1, indices_dict=indices_dict)

        return coord_dicts

def get_coordinates(traj_and_report, indices_dict):
    trajfile, report = traj_and_report
    indices = np.concatenate([i[0] for i in indices_dict.values()])
    accepted_steps = hl.accepted_pele_steps(report)

    traj = hl.load_trajectory(trajfile, indices)
    coords = traj.xyz *10  # coord units from nm to A
    #coords = coords[accepted_steps] # duplicate rows when a step is rejected

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
            frag_regex = ".*(?P<frag>frag\d+$)"
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
        frag_regex = ".*(?P<frag>frag\d+$)"
        frag = re.match(frag_regex, d)['frag']
        features = fragment_features[frag]
        return features
