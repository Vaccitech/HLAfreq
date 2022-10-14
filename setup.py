import setuptools

setuptools.setup(
    name="HLAfreq",
    version="0.0.1dev0",
    url="https://github.com/Vaccitech/HLAfreq",
    project_urls={
        'Tracker': "https://github.com/Vaccitech/HLAfreq/issues"
    },
    author="David Wells",
    author_email="david.wells@vaccitech.co.uk",
    description="Download and combine HLA frequency data from multiple studies",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=[
        'bs4',
        'requests',
        'pandas',
        'numpy',
        'matplotlib',
        'scipy',
    ],
    classifiers=[
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
    ],
    include_package_data=True
)