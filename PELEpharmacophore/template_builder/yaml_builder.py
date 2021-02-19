from string import Template
from PELEpharmacophore.template_builder import base_class as bc

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


class YamlBuilder(bc.TemplateBuilder):
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
