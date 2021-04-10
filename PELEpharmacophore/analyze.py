import os
import re
from argparse import ArgumentParser
from PELEpharmacophore.errors import custom_errors as ce
import PELEpharmacophore.analysis.grid_analyzer as ga
import PELEpharmacophore.analysis.meanshift_analyzer as sa
import PELEpharmacophore.yaml_parser as yp
import PELEpharmacophore.valid_flags as vf
import PELEpharmacophore.data.fragment_features as ff


def parse_args():
    '''
    Command line parser
    '''
    parser = ArgumentParser(description='Run PELEpharmacophore analysis')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args()
    return args.input_file

def run_PELEpharmacophore(analyzer, dir, chain, resname, resnum, center, radius, features, ncpus):
    analyzer.set_dir(dir)
    analyzer.set_ligand(chain, resname, resnum)
    analyzer.set_features(features)
    if center:
        analyzer.set_grid(center, radius)
    analyzer.run(ncpus)


def PELEpharmacophore_ligand(analyzer_class, dir, chain, resname, resnum, center, radius, features, ncpus, filt=2):
    analyzer = analyzer_class(dir)
    analyzer = run_PELEpharmacophore(dir, chain, resname, resnum, center, radius, features, ncpus)
    analyzer.set_frequency_filter(filt) # add filter as arg
    analyzer.save_pharmacophores()
    return

def PELEpharmacophore_fragments(analyzer_class, dir, center, radius, ncpus=20, fragment_features=ff.fragment_features, filt=0):
    analyzer = analyzer_class(dir)
    subdirs = [os.path.join(dir, subdir) for subdir in os.listdir(dir)]

    for frag_dir in subdirs:
        frag_regex = ".*(?P<frag>frag\d+$)"
        frag = re.match(frag_regex, frag_dir)['frag']
        features = fragment_features[frag]
        analyzer = run_PELEpharmacophore(analyzer, frag_dir, "L", "FRA", 900, center, radius, features, ncpus)


    analyzer_frags.set_frequency_filter(filt)
    analyzer_frags.save_pharmacophores("Pharmacophores_frags")
    return

analysis_types = {
                 "grid": ga.GridAnalyzer,
                 "meanshift": ma.MeanshiftAnalyzer
}

def main(input_yaml):
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS)
    try:
        yaml_obj.read()
    except AttributeError:
        raise ce.WrongYamlFile(f"Input file: {input_yaml} does not look like a correct yml file")

    analysis_class = analysis_types[yaml_obj.analysis_type]

    if yaml_obj.ligand:
        PELEpharmacophore_ligand(analysis_class, yaml_obj.dir, yaml_obj.chain, yaml_obj.resname, yaml_obj.resnum, yaml_obj.grid_center, yaml_obj.grid_radius, yaml_obj.features)
    else:
        PELEpharmacophore_fragments(analysis_class, yaml_obj.dir, yaml_obj.grid_center, yaml_obj.grid_radius)


if __name__ == "__main__":
    args = parse_args()
    main(args)
