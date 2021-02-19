import os
import re
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

    def __init__(self, target, ligand, center, final_dir="systems"):
        self.target = target
        self.ligand = ligand
        self.center = center
        self.final_dir = final_dir

    def pdbconvert(self, file, in_format="pdb", out_format="mae", outdir="."):
        input_path = os.path.abspath(file)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
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
        glide_job, ext = os.path.splitext(self.glide_input)
        self.glide_output = f"{glide_job}_pv.maegz"
        os.system(command)

    def rename_files(self, indir = "docking_results"):
        self.docking_dir = indir
        filelist =[os.path.join(self.docking_dir, f) for f in os.listdir(self.docking_dir)] #poner en función aparte
        pattern = "^TITLE {5}(.*)"
        for f in filelist:
            if f.endswith("1.pdb"):
                os.rename(f, self.target)
            else:
                with open(f, 'r') as fh:
                    for line in fh:
                         match = re.search(pattern, line)
                         if match:
                             pdb_name = match.group(1)
                             break
                os.rename(f, f"{os.path.join(indir, pdb_name)}.pdb")

    def create_systems(self, ligchain="L", ligname="FRA"):
        if not os.path.isdir(self.final_dir):
            os.mkdir(self.final_dir)

        filelist =[os.path.join(self.docking_dir, f) for f in os.listdir(self.docking_dir)]
        targetname, ext = os.path.splitext(self.target)

        for i, f in enumerate(filelist):
            fragname, ext = os.path.splitext(os.path.basename(f))
            system = open(f"{self.final_dir}/{targetname}_{fragname}.pdb", "w")
            with open(self.target, "r") as targetfile:
                [system.write(line) for line in targetfile if not line.startswith("END")]
            with open(f, "r") as fragfile:
                for line in fragfile:
                    if line.startswith(("HETATM", "CONECT")):
                        line = line.replace("UNK   ", f"{ligname} L ")
                        system.write(line)
                system.write("ENDMDL\nEND\n")
            system.close()

    def run(self):
        self.pdbconvert(self.target)
        while not os.path.exists(self.convert_output):
            time.sleep(10)
        self.generate_glide_grids()
        while not os.path.exists(self.gridfile):
            time.sleep(10)
        self.generate_glide_input()
        self.glide()
        while not os.path.exists(self.glide_output):
            time.sleep(10)
        self.pdbconvert(self.glide_output, in_format="mae", out_format="pdb", outdir="docking_results") #!!!! hay que indicar input file
        self.rename_files()
        self.create_systems()


if __name__ == "__main__":
    args = parse_args()
    docking = GlideDocking(args.target, args.ligand, args.center)
    docking.run()
