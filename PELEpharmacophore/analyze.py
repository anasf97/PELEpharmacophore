import os
import re
from argparse import ArgumentParser
from PELEpharmacophore.errors import custom_errors as ce
import PELEpharmacophore.analysis.grid_analyzer as ga
import PELEpharmacophore.analysis.meanshift_analyzer as ma
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

def run_PELEpharmacophore(analyzer, indir, chain, resname, resnum, center, radius, features, ncpus):
    analyzer.set_dir(indir)
    analyzer.set_ligand(chain, resname, resnum)
    analyzer.set_features(features)
    analyzer.run(ncpus)

def PELEpharmacophore_ligand(analyzer_class, indir, chain, resname, resnum, center, radius, features, ncpus, outdir, filt=2):
    analyzer = analyzer_class(indir)
    analyzer = run_PELEpharmacophore(indir, chain, resname, resnum, center, radius, features, ncpus)
    analyzer.set_frequency_filter(filt) # add filter as arg
    analyzer.save_pharmacophores(outdir)

def PELEpharmacophore_fragments(analyzer_class, indir, center, radius, outdir, ncpus=23, fragment_features=ff.fragment_features, filt=1):
    import datetime

    start = datetime.datetime.now()
    t = start.strftime("%H:%M:%S")
    print("Starting time:", t)

    analyzer = analyzer_class(indir)
    if isinstance(analyzer, ga.GridAnalyzer):
        analyzer.set_grid(center, radius)

    analyzer.run(ncpus)
    analyzer.set_frequency_filter(filt)
    analyzer.save_pharmacophores(outdir)

    finish = datetime.datetime.now()
    t = finish.strftime("%H:%M:%S")
    print("Finishing time:", t)

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
        PELEpharmacophore_ligand(analysis_class, yaml_obj.dir, yaml_obj.chain, yaml_obj.resname, yaml_obj.resnum, yaml_obj.grid_center, yaml_obj.grid_radius, yaml_obj.features, yaml_obj.outdir)
    else:
        PELEpharmacophore_fragments(analysis_class, yaml_obj.dir, yaml_obj.grid_center, yaml_obj.grid_radius, yaml_obj.outdir)


if __name__ == "__main__":
    args = parse_args()
    main(args)
