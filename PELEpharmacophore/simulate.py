import argparse
from PELEpharmacophore.simulation import docking as dk
from PELEpharmacophore.simulation import launch_file_creator as lf
from PELEpharmacophore.simulation import simulation_runner as sr

def parse_args():
    parser = argparse.ArgumentParser(description="Glide docking for a pdb target file and a mae ligand file.")
    parser.add_argument("target", help="PDB apo target file")
    parser.add_argument("ligand", help="Mae ligand file")
    parser.add_argument("center", help="Center of the binding site")
    args = parser.parse_args()
    return args

def main(target, ligand, center):
    docking = dk.GlideDocking(target, ligand, center)
    launch_file_creator = lf.LaunchFileCreator(docking.final_dir)
    sr.SimulationRunner(launch_file_creator.slurm_outdir)

if __name__ == '__main__':
    args = parse_args(arg.target, args.ligand, args.center)
    main()
