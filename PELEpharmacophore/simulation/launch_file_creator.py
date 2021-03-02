import os
import argparse
from PELEpharmacophore.template_builder import yaml_builder as yb
from PELEpharmacophore.template_builder import slurm_builder as sb

def parse_args():
    parser = argparse.ArgumentParser(description="Create slurm and yaml file for all systems in a directory.")
    parser.add_argument("indir", help="Input directory.")
    args = parser.parse_args()
    return args

class LaunchFileCreator(object):
    """docstring for LaunchFileCreator."""

    def __init__(self, indir, ncpus = 24, ligchain = "L", ligname = "FRA", simulation_folder="PELEpharmacophore_simulations"):
        self.filelist = [os.path.join(indir, f) for f in os.listdir(indir)]
        self.ligchain = ligchain
        self.ligname = ligname
        self.ncpus = ncpus
        self.simulation_folder = simulation_folder

        self.create_launch_files()

    def create_launch_files(self, yaml_outdir="yaml_files", slurm_outdir="slurm_files"):
        self.slurm_outdir = slurm_outdir
        for f in self.filelist:
            system_pdb = os.path.basename(f)
            system, ext = os.path.splitext(system_pdb)

            working_folder = os.path.join(self.simulation_folder, system)
            
            yaml_args = self.yaml_args(f, self.ligchain, self.ligname, working_folder)
            yaml_name = system_pdb.replace(".pdb", ".yml")
            yaml_path = os.path.join(yaml_outdir, yaml_name)
            yb.YamlBuilder(yaml_args, yaml_name, yaml_outdir)

            slurm_args = self.slurm_args(working_folder, yaml_path, self.ncpus)
            slurm_name = system_pdb.replace(".pdb", ".sl")
            sb.SlurmBuilder(slurm_args, slurm_name, slurm_outdir)

    def yaml_args(self, system, ligchain, ligname, working_folder):
        return {"system" : system,
                "chain" : ligchain,
                "resname" : ligname,
                "working_folder" : working_folder
                }

    def slurm_args(self, system, system_yml, ncpus):
        return {"system" : system,
                "system_yml" : system_yml,
                "ncpus" : ncpus,
                }

if __name__ == "__main__":
    args = parse_args()
    LaunchFileCreator(args.indir)
