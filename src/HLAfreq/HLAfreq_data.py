"""
Data loaders
"""

import pkg_resources
import pandas as pd


def load_countries():
    """Load dataframe of countries and regions available on [allele frequencies.net](http://www.allelefrequencies.net/).

    Returns:
        pd.DataFrame: Three columns: `Country`, `Region`, and `largeRegion`.
        `largeRegion` is based on `Region` with with fewer categories to improve
        plotting.
    """
    stream = pkg_resources.resource_stream(__name__, "data/countries.csv")
    return pd.read_csv(stream, encoding="latin-1")


def load_HLA1supertypes_Sidney2008():
    """Load HLA alleles and their supertype as defined in Sidney2008

    Returns:
        pd.DataFrame: DataFrame of alleles and their supertype
    """
    stream = pkg_resources.resource_stream(
        __name__, "data/HLA1supertypes_Sidney2008.csv"
    )
    return pd.read_csv(stream, encoding="latin-1")
