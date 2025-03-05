import setuptools
import warnings

requirements = ["besmarts", 'numpy']

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='besmarts-mechanics',
    version='0.1.1',
    description="BESMARTS molecular mechanics force field fitting",
    license="MIT",
    author="Trevor Gokey",
    author_email='tgokey@uci.edu',
    url='https://github.com/trevorgokey/besmarts',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=[
        'besmarts',
        'besmarts.mechanics'
    ],
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

warnings.warn("""
    WARNING!!! The default behavior of loading SMIRNOFF force fields uses
    the sqm program from Ambertools to determine partial charges. Ambertools
    is not available on PyPI at the time of this writing. Please install
    Ambertools using one of the supported methods and ensure that sqm is in
    your PATH. Binary packages are available through the conda package manager
    via
    
    conda install -c conda-forge ambertools
    
    """,
    RuntimeWarning
)
