import os
from string import Template

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
