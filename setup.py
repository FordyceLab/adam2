from setuptools import setup

setup(
    name='adam2',
    version='0.1',
    description='Hybridization oligo design utility',
    url='https://github.com/FordyceLab/adam2',
    author='Tyler Shimko',
    license='MIT',
    packages=['adam2'],
    install_requires=[
        'biopython',
        'distance',
        'numpy',
        'primer3-py',
        'tqdm',
        'pyaml',
    ],
    scripts=['bin/adam2']
)