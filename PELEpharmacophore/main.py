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

def run_PELEpharmacophore(dir, chain, resname, resnum, center, features):
    target = tr.Target(dir)
    target.set_ligand(chain, resname, resnum)
    target.set_features(features)
    target.set_grid(center)
    with Pool(1) as p: # add n_workers as arg
        dicts = p.map(target.analyze_trajectory, target.filelist)
        p.close()
        p.join()
    for d in dicts:
        target.merge_voxel_dicts(d)
    target.get_frequencies()
    target.set_frequency_filter(2) # add filter as arg
    target.save_pharmacophores()

def analyze_fragment_simulations(dir, center):
    subdirs = os.listdir(dir)
    for i, frag in enumerate(fragments):
        run_PELEpharmacophore(subdir[i], "L", FR, 900, )
    pass

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS)
    #try:
    yaml_obj.read()
    #except AttributeError:
        #raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")
    run_PELEpharmacophore(yaml_obj.dir, yaml_obj.chain, yaml_obj.resname, yaml_obj.resnum, yaml_obj.grid_center, yaml_obj.features)


if __name__ == "__main__":
    args = parse_args()
    main(args)
