from dataclasses import dataclass
from difflib import SequenceMatcher
import os
import yaml
import warnings


@dataclass
class YamlParser(object):

    yamlfile: str
    valid_flags: dict

    def read(self) -> None:
        self.data = self._parse_yaml()
        self._check()
        self._parse()

    def _parse_yaml(self) -> dict:
        # Retrieve raw info from yaml
        with open(self.yamlfile, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise(exc)
        return data

    def _check(self) -> None:
        # Check if valids in yaml file are valids
        for key in self.data.keys():
            if key not in self.valid_flags.values():
                raise KeyError(self._recommend(key))

    def _recommend(self, key):
        most_similar_flag = None
        for valid_key in self.valid_flags.values():
            flag = Most_Similar_Flag(valid_key)
            flag.calculate_distance(key)
            if not most_similar_flag:
                most_similar_flag = flag
            else:
                 if flag.distance > most_similar_flag.distance:
                     most_similar_flag = flag
        exception_raised = f"Incorrect flag {key}. Did you mean {most_similar_flag.name}?"
        return exception_raised

    def _parse(self) -> None:
        # Parse fields in yaml file and set defaults
        valid_flags = self.valid_flags
        data = self.data
        self.dir = data.get(valid_flags["dir"], "")
        self.dir = os.path.abspath(self.dir) if self.dir else ""
        self.chain = data.get(valid_flags["chain"], "")
        self.resname = data.get(valid_flags["resname"], "")
        self.resnum = data.get(valid_flags["resnum"], "")
        self.features = data.get(valid_flags["features"], None)
        self.grid_center = data.get(valid_flags["grid_center"], None)
        self.grid_radius = data.get(valid_flags["grid_radius"], None)
        self.ligand = data.get(valid_flags["ligand"], None)
        #self.verbose = data.get(valid_flags["verbose"], None)

@dataclass
class Most_Similar_Flag():

    name: str

    def calculate_distance(self, key):
        self.distance = SequenceMatcher(None, self.name, key).ratio()
