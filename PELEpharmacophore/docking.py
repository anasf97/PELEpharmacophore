import os
import time
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Glide docking for a pdb target file and a mae ligand file.")
    parser.add_argument("target", help="PDB apo target file")
    parser.add_argument("ligand", help="Mae ligand file")
    parser.add_argument("center", help="Center of the binding site")
    args = parser.parse_args()
    return args

class GlideDocking:

    def __init__(self, target, ligand, center):
        self.target = target
        self.ligand = ligand
        self.center = center

    def pdbconvert(self, in_format="pdb", out_format="mae", outdir="."):
        input_path = os.path.abspath(self.target)
        self.convert_output = os.path.basename(self.target.replace(f".{in_format}", f".{out_format}"))
        self.convert_output = os.path.join(outdir, self.convert_output)
        schrodinger_path ="$SCHRODINGER/utilities/pdbconvert"
        command = f"{schrodinger_path} -i{in_format} {input_path} -o{out_format} {self.convert_output}"
        os.system(command)

    def generate_glide_grids(self):
        schrodinger_path = "$SCHRODINGER/utilities/generate_glide_grids"
        command = f"{schrodinger_path}  -rec_file {self.convert_output} -cent_coor '{self.center}'"
        self.gridfile = "generate-grids-gridgen.zip"
        os.system(command)

    def generate_glide_input(self, precision ="SP", poses="1", glide_input ="input_glide.in"):
        lines = [f"GRIDFILE {self.gridfile}",
                 f"LIGANDFILE {self.ligand}",
                 f"POSES_PER_LIG {poses}",
                 f"PRECISION {precision}"]
        self.glide_input = glide_input
        with open(glide_input, "w") as f:
            f.writelines(f"{l}\n" for l in lines)

    def glide(self):
        schrodinger_path = f"$SCHRODINGER/glide"
        command = f"{schrodinger_path} {self.glide_input}"
        os.system(command)

    def run(self):
        self.pdbconvert()
        while not os.path.exist(self.convert_output):
            time.sleep(10)
        self.generate_glide_grids()
        while not os.path.exist(self.gridfile):
            time.sleep(10)
        self.generate_glide_input()
        self.glide()
        self.pdbconvert(in_format="mae", out_format="pdb", outdir="docking_results")

if __name__ == "__main__":
    args = parse_args()
    docking = GlideDocking(args.target, args.ligand, args.center)
    docking.run()
