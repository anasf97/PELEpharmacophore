from argparse import ArgumentParser
from multiprocessing import Pool
import PELEpharmacophore.target as tr
import PELEpharmacophore.yaml_parser as yp
import PELEpharmacophore.valid_flags as vf
import PELEpharmacophore.frag_features as ff

fragments = [

]
def parse_args(args=[]):
    '''
    Command line parser
    '''
    parser = ArgumentParser(description='Run PELE Platform')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args.input_file

def run_PELEpharmacophore(dir, chain, resname, resnum, center, radius, features):
    target = tr.Target(dir)
    target.set_ligand(chain, resname, resnum)
    target.set_features(features)
    target.set_grid(center, radius)
    with Pool(1) as p: # add n_workers as arg
        grids = p.map(target.analyze_trajectory, target.filelist)
        p.close()
        p.join()
    for g in grids:
        target.merge_grids(g)
    return target

def PELEpharmacophore_ligand(dir, chain, resname, resnum, center, radius, features):
    target = run_PELEpharmacophore
    target.set_frequency_filter(2) # add filter as arg
    target.save_pharmacophores()
    return target

def PELEpharmacophore_fragments(dir, center, radius):
    target_frags = Target()
    target_frags.set_grid(center, radius)
    subdirs = os.listdir(dir)
    for i, features in enumerate(fragments):
        target = run_PELEpharmacophore(subdir[i], "L", "FRA", 900, center, radius, features)
        target_frags.merge_grids(target.grid)
    target_frags.set_frequency_filter(2)
    target_frags.save_pharmacophores()
    return target_frags

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS)
    #try:
    yaml_obj.read()
    #except AttributeError:
        #raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")
    run_PELEpharmacophore(yaml_obj.dir, yaml_obj.chain, yaml_obj.resname, yaml_obj.resnum, yaml_obj.grid_center, yaml_obj.grid_radius, yaml_obj.features)


if __name__ == "__main__":
    args = parse_args()
    main(args)
