import os

def parse_args():
    parser = argparse.ArgumentParser(description="Launch simlations from all slurm files in a directory.")
    parser.add_argument("indir", help="Input directory.")
    args = parser.parse_args()
    return args

class SimulationRunner(object):
    """docstring for SimulationRunner."""

    def __init__(self, indir):
        self.filelist = [os.path.join(indir, f) for f in os.listdir(indir)]
        self.run()

    def sbatch(self, run_file):
        command = f"sbatch {run_file}"
        os.system(command)

    def run(self):
        for f in self.filelist:
            self.sbatch(f)

if __name__ == '__main__':
    args = parse_args()
    SimulationRunner(args.indir)
