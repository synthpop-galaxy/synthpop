import os

from setuptools import setup, find_packages, Command

setup(
    name='synthpop',
    version='0.1.0',
    description='SynthPop is a modular framework to generate synthetic population models',
    author='J. Klüter, S. Johnson, M.J. Huston, A. Aronica, M. Penny',
    license='GPLv3',
    packages=find_packages(),

    include_package_data=True,
    install_requires=[
        'numpy>=1.23.2',
        'scipy>=1.9.0',
        'pandas>=1.4.0',
        'requests>=2.28.1',
        'tables>=3.7.0',
        'tqdm~=4.64.1'
    ],
    extras_require={
        'optional': ['astropy>=4.3.1', 'dustmaps>=1.0.10', 'astroquery~=0.4.6']
    },
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown'
    cmdclass = {
    'uninstall': UninstallCommand,
    },
)

class UninstallCommand(Command):
    description = 'Uninstall the package and remove any generated symbolic links.'
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Remove any generated symbolic links."""
        # Remove file1.txt
        if os.path.islink("modules"):
            os.unlink(symlink)        # Remove file2.txt
        if os.path.islink("models"):
            os.unlink("models")
        if os.path.islink("config_files"):
            os.unlink("config_files")
        if os.path.islink("constants.py"):
            os.unlink('constants.py')
