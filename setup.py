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
    version='1.5.3',
    description='Python library containing tools for dealing with conversions between particle energy and rigidity',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Chris S. W. Davis',
    author_email='ChrisSWDavis@gmail.com',
    license='MIT',
    url='https://github.com/ChrisSWDavis/ParticleRigidityCalculationTools',
    install_requires=install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
