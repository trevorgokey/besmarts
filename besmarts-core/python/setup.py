import setuptools

requirements = []

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='besmarts',
    version='0.1.0',
    description="Binary Encoded SMARTS",
    license="MIT",
    author="Trevor Gokey",
    author_email='tgokey@uci.edu',
    url='https://github.com/trevorgokey/besmarts',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=[
        'besmarts',
        'besmarts.assign',
        'besmarts.cluster',
        'besmarts.codecs',
        'besmarts.core',
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
    extras_require={
        "rdkit": ['besmarts-rdkit'],
        "scipy": ['besmarts-scipy'],
        "openmm": ['besmarts-openmm'],
    }
)
