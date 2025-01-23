from setuptools import find_packages, setup
import os

# get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='ParticleRigidityCalculationTools',
    py_modules=["ParticleRigidityCalculationTools"],
    version='1.5.17',
    description='Python library containing tools for dealing with conversions between particle energy and rigidity',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Space Environment and Protection Group, University of Surrey',
    url='https://github.com/ssc-maire/ParticleRigidityCalculationTools',
    keywords = 'space physics earth geomagnetic rigidity magnetocosmics',
    license='CC BY-NC-SA 4.0',
    install_requires=['numpy>=1.23.1',
                      'pandas>=1.4.3',
                      'setuptools>=45.2.0'],
    #install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
