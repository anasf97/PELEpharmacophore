import os
from string import Template
from PELEpharmacophore.template_builder import base_class as bc

PYTHONPATH = os.environ.get("PYTHONPATH")
SLURM_TEMPLATE = f"{PYTHONPATH}/PELEpharmacophore/templates/slurm_template.sl"

slurm_default_values = {"system": None,
                        "system_yml":None,
                        "ncpus": 24,
                        "pele_path": "/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/",
                        "schrodinger_path": "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC",
                        "conda_path": "/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin",
                        }

class SlurmBuilder(bc.TemplateBuilder):

    default_values = slurm_default_values
    template = SLURM_TEMPLATE

    def __init__(self, keywords, filename= "input.sl", outdir="."):
        super().__init__(keywords, filename, outdir)

    def fill_in(self):
        self.default_values.update(self.keywords)
        with open(self.template, 'r') as infile:
            confile_data = infile.read()

        confile_template = Template(confile_data)

        self.lines = confile_template.safe_substitute(self.default_values)
