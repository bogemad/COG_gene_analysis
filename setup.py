from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "COG_gene_analysis",
    version = "0.1",
    author = "Daniel Bogema",
    author_email = "daniel.bogema@dpi.nsw.gov.au",
    description = ("Scripts to create graphs showing numbers of COG genes and identify unique COG proteins in comparative genomics"),
    license = "GPL-3.0",
    keywords = "genomics gene function",
    url = "https://github.com/bogemad/COG_gene_analysis",
    py_modules=['COG_gene_analysis/COGs', 'COG_gene_analysis/unique_proteins', 'COG_gene_analysis/venn_diagram_input'],
    scripts=['bin/generate_figures_from_orthofinder_output'],
    long_description=read('README.md'),
)