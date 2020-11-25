from argparse import ArgumentParser
from multiprocessing import Pool
import PELEpharmacophore.target as tr
import PELEpharmacophore.yaml_parser as yp

VALID_FLAGS = { "dir": "dir",
                "chain":"chain",
                "resname": "resname",
                "resnum": "resnum",
                "HBD": "HBD",
                "HBA": "HBA",
                "ARO": "ARO",
                "ALI": "ALI",
                "NEG": "NEG",
                "POS": "POS",
                "features":"features"}

def parse_args(args=[]):
    '''
    Command line parser
    '''
    parser = ArgumentParser(description='Run PELE Platform')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args.input_file

def run_PELEpharmacophore(args):
    target = tr.Target(args.dir)
    target.set_ligand(args.chain, args.resname, args.resnum)
    #features = {"HBD":args.HBD, "HBA":args.HBA, "ARO":args.ARO, "ALI":args.ALI, "NEG":args.NEG, "POS":args.POS}
    features = args.features
    target.set_features(features)
    target.set_grid((2.173, 15.561, 28.257), 7)
    with Pool(1) as p: # add n_workers as arg
        dicts = p.map(target.analyze_trajectory, target.filelist)
        p.close()
        p.join()
    for d in dicts:
        target.merge_voxel_dicts(d)
    target.get_frequencies()
    target.set_frequency_filter(2) # add filter as arg
    target.save_pharmacophores()

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, VALID_FLAGS)
    #try:
    yaml_obj.read()
    #except AttributeError:
        #raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")
    pass
    run_PELEpharmacophore(yaml_obj)


if __name__ == "__main__":
    input_file = parse_args()
    main(input_file)
