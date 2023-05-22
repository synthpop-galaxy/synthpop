from setuptools import setup, find_packages

setup(
    name='synthpop',
    version='0.2.0',
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
        'tqdm~=4.64.1',
        'pydantic~=1.10.7'
        ],
    extras_require={
        'optional': ['astropy>=4.3.1', 'dustmaps>=1.0.10', 'astroquery~=0.4', 'matplotlib~=3.6.2']
        },
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown'
    )
