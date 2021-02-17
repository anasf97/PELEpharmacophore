import os
from string import Template


yaml_default_values = {"system": None,
              "chain": None,
              "resname": None,
              "induced_fit_exhaustive": None,
              "seed": 12345,
              "usesrun": None,
              "working_folder": None,
              "cpus": None,
              "pele_license": "/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license",
              "pele_exec": "/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi",
            }

YAML_TEMPLATE = "$flag: $value\n"


class TemplateBuilder(object):

    def __init__(self, keywords, filename, outdir):
        self.outdir = outdir
        self.filename = filename
        self.keywords = {k: v for k, v in keywords.items() if v is not None}

        self.fill_in()
        self.save()

    def save(self):
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        outfile = os.path.join(self.outdir, self.filename)
        with open(outfile, 'w') as of:
            of.writelines(self.lines)


class YamlBuilder(TemplateBuilder):
    default_values = yaml_default_values
    template = YAML_TEMPLATE

    def __init__(self, keywords, filename="input.yaml", outdir="."):
        super().__init__(keywords, filename, outdir)

    def fill_in(self):
        yaml_template = Template(self.template)
        self.default_values.update(self.keywords)
        self.lines = []
        for k, v in self.default_values.items():
            if v is not None:
                 self.lines.append(yaml_template.safe_substitute(flag = k, value = v))


slurm_default_values = {"system": None,
                          "ncpus": 24,
                          "pele_path": "/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/",
                          "schrodinger_path": "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC",
                          "conda_path": "/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin",
                          }

class SlurmBuilder(TemplateBuilder):
    default_values = slurm_default_values
    template = "slurm_template.sl"

    def __init__(self, keywords, filename= "input.sl", outdir="."):
        super().__init__(keywords, filename, outdir)


    def fill_in(self):
        self.default_values.update(self.keywords)
        with open(self.template, 'r') as infile:
            confile_data = infile.read()

        confile_template = Template(confile_data)

        self.lines = confile_template.safe_substitute(self.keywords)


if __name__ == "__main__":
    args = {"system": 'systems/apo_1a9u_frag0.pdb',
            "chain": 'L',
            "induced_fit_exhaustive": "true",
            "resname": 'FRA',
            }
    args2 = {"system": "apo_1a9u.pdb",
              "ncpus": 24,
              "pele_path": "/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/",
              "schrodinger_path": "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC",
              "conda_path": "/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin",
              }
    SlurmBuilder(args2)
