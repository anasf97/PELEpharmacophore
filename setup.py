import setuptools


setuptools.setup(name='PELEpharmacophore',
      version="1.0.0",
      url = "https://github.com/anasf97/PELEpharmacophore",
      description='Generate pharmacophore models.',
      author='Ana Sánchez',
      author_email='ana.sanchez01@estudiant.upf.edu',
      packages=setuptools.find_packages(),
      install_requires=["biopython", "scipy", "dataclasses", "sklearn", "pytest", "flake8", "mdtraj"],
      classifiers=[
       "Programming Language :: Python :: 3",
       "License :: OSI Approved :: MIT License",
       "Operating System :: OS Independent",
       ],
     )
