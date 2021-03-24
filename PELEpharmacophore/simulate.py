from argparse import ArgumentParser
from PELEpharmacophore.errors import custom_errors as ce
from PELEpharmacophore.simulation import docking as dk
from PELEpharmacophore.simulation import launch_file_creator as lf
from PELEpharmacophore.simulation import simulation_runner as sr
import PELEpharmacophore.yaml_parser as yp
import PELEpharmacophore.valid_flags as vf

def parse_args():
    """
    Command line parser
    """
    parser = ArgumentParser(description='Run PELEpharmacophore simulations')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args()
    return args.input_file

def simulate(target, ligand, center):
    #docking = dk.GlideDocking(target, ligand, center)
    #docking.run()
    #launch_file_creator = lf.LaunchFileCreator(docking.final_dir)
    launch_file_creator = lf.LaunchFileCreator("subset/systems")
    sr.SimulationRunner(launch_file_creator.slurm_outdir)

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS)
    try:
        yaml_obj.read()
    except AttributeError:
        raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")
    simulate(yaml_obj.target, yaml_obj.ligand, yaml_obj.grid_center)

if __name__ == '__main__':
    args = parse_args()
    main(args)
