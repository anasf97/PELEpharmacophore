import os
import time
import numpy as np
import shutil as sh
import PELEpharmacophore.helpers as hl


PYTHONPATH = os.environ.get("PYTHONPATH")
FEATURE_DEF =  f"{PYTHONPATH}/PELEpharmacophore/templates/pharma_feature.ini"

EQ_FEATURES = {
                "HBA" : "A",
                "HBD" : "D",
                "ARO" : "R",
                "ALI" : "H",
                "NEG" : "N",
                "POS" : "P"
}

class PharmacophoreWriter(object):
    """docstring for PharmacophoreWriter."""

    def __init__(self, name, feature_dict, coords, outdir):
        self.name = os.path.basename(name)
        updated_dict = {EQ_FEATURES[feature] : value for feature, value in feature_dict.items()}
        self.feature_dict = updated_dict
        self.coords = coords
        self.outdir = outdir

        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

        self.run()


    def write_files(self, tol=2):
        xyz_lines, dxyz_lines, rules_lines, mask_lines = ([] for i in range(4))

        for feature, results in sorted(self.feature_dict.items()):
            for i, site in enumerate(results):
                x, y, z = site.center
                i += 1
                xyz = f"{i} {feature} {x:>7.4f} {y:>7.4f} {z:>7.4f}\n"
                dxyz = f"{i} {feature} {tol}\n"
                rules = f"{i} {feature}"
                mask = f"{i} {feature} 0"

                xyz_lines.append(xyz)
                dxyz_lines.append(dxyz)
                rules_lines.append(rules)
                mask_lines.append(mask)

        all_lines = [xyz_lines, dxyz_lines, rules_lines, mask_lines]
        extensions = ["xyz", "dxyz", "rules", "mask"]

        filenames = [hl.custom_path(".", self.name, ext) for ext in extensions]
        print(filenames)

        for lines, filename in zip(all_lines, filenames):
            self.save(lines, filename)


    def write_def(self):
        filename = f"{self.name}.def"
        sh.copyfile(FEATURE_DEF, filename)


    def write_shell(self):
        shell = os.path.join(self.outdir, "shell.pdb")
        with open(shell, 'w') as f:
            shell_coords = np.unique(self.coords, axis=0)
            for i, point in enumerate(shell_coords):
                print(point)
                f.write(hl.format_line_pdb(point, "C", resnum=i))

        hl.pdbconvert(shell, outdir=self.outdir)


    def save(self, lines, filename):
        outfile = os.path.join(self.outdir, filename)
        with open(outfile, 'w') as of:
            of.writelines(lines)


    def generate_pharmacophore(self):
        schrodinger_path ="$SCHRODINGER/utilities/phase_hypo_util"
        command = f"{schrodinger_path} {self.name} pack"
        os.system(command)


    def generate_xvols(self):
        schrodinger_path ="$SCHRODINGER/utilities/create_xvolShell"
        mae = os.path.join(self.outdir, "shell.mae")
        command = f"{schrodinger_path} {self.name} -ref {mae} -buff 2"
        os.system(command)


    def run(self):
        output = os.path.join(self.outdir, f"{self.name}.phypo")
        self.write_files()
        self.write_def()
        self.write_shell()
        self.generate_pharmacophore()
        while not os.path.exists(output):
            time.sleep(1)
        generate_xvols()
