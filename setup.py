from setuptools import setup, find_packages

setup(
    name="HLAfreq",
    version="0.0.1dev3",
    url="https://github.com/Vaccitech/HLAfreq",
    project_urls={
        'Tracker': "https://github.com/Vaccitech/HLAfreq/issues"
    },
    author="David Wells",
    author_email="david.wells@vaccitech.co.uk",
    description="Download and combine HLA frequency data from multiple studies",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=[
        'bs4',
        'requests',
        'pandas',
        'numpy',
        'matplotlib',
        'scipy',
        'pymc'
    ],
    classifiers=[
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
    ],
    package_dir={'':'src'},
    packages=find_packages(
        where='src'
    ),
    include_package_data=True,
)