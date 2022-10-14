"""
Data loaders
"""
import pkg_resources
import pandas as pd

def load_countries():
    stream = pkg_resources.resource_stream(__name__, 'data/countries.csv')
    return pd.read_csv(stream, encoding='latin-1')

def load_HLA1supertypes_Sidney2008():
    stream = pkg_resources.resource_stream(__name__, 'data/HLA1supertypes_Sidney2008.csv')
    return pd.read_csv(stream, encoding='latin-1')
