from argparse import ArgumentParser
from PELEpharmacophore.errors import custom_errors as ce
import PELEpharmacophore.analysis.simulation_analyzer as sa
import PELEpharmacophore.yaml_parser as yp
import PELEpharmacophore.valid_flags as vf
import PELEpharmacophore.data.fragment_features as ff


def parse_args(args):
    '''
    Command line parser
    '''
    parser = ArgumentParser(description='Run PELEpharmacophore analysis')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args.input_file

def run_PELEpharmacophore(dir, chain, resname, resnum, center, radius, features, ncpus):
    analyzer = sa.SimulationAnalyzer(dir)
    analyzer.set_ligand(chain, resname, resnum)
    analyzer.set_features(features)
    analyzer.set_grid(center, radius)
    analyzer.run(ncpus)
    return analyzer

def PELEpharmacophore_ligand(dir, chain, resname, resnum, center, radius, features, ncpus, filter=2):
    analyzer = run_PELEpharmacophore(dir, chain, resname, resnum, center, radius, features, ncpus)
    analyzer.set_frequency_filter(filter) # add filter as arg
    analyzer.save_pharmacophores()
    return analyzer

def PELEpharmacophore_fragments(dir, center, radius, ncpus=20, fragment_features=ff.fragment_features, filter=2):
    analyzer_frags = SimulationAnalyzer()
    analyzer_frags.set_grid(center, radius)
    subdirs = os.listdir(dir)
    for frag, features in fragment_features.items():
        frag_dir = [s for s in subdirs if s.endswith(frag)][0]
        analyzer = run_PELEpharmacophore(frag_dir, "L", "FRA", 900, center, radius, ncpus, features)
        analyzer_frags.merge_grids(analyzer_frags.grid, analyzer.grid)
    analyzer_frags.set_frequency_filter(filter)
    analyzer_frags.save_pharmacophores()
    return analyzer_frags

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS)
    try:
        yaml_obj.read()
    except AttributeError:
        raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")
    if yaml_obj.ligand:
        PELEpharmacophore_ligand(yaml_obj.dir, yaml_obj.chain, yaml_obj.resname, yaml_obj.resnum, yaml_obj.grid_center, yaml_obj.grid_radius, yaml_obj.features)
    else:
        PELEpharmacophore_fragments(yaml_obj.dir, yaml_obj.grid_center, yaml_obj.grid_radius)


if __name__ == "__main__":
    args = parse_args()
    main(args)
