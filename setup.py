#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:58:20 2020

@author: rasmus
"""

import setuptools
import pandas as pd
import os

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'


# First we assemble the correlation matrix
picklepath = 'data/buildpickles/'
pfiles = os.listdir(picklepath)
dfs = []
for pickle_file in pfiles:
    dfs.append(pd.read_pickle(picklepath + pickle_file))
dfs = pd.concat(dfs).sort_index()
dfs.to_pickle('data/pickles/correlations.p')


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TFenrich",
    version="0.01",
    author="Rasmus Magnusson",
    author_email="rasma774@gmail.com",
    description="A package to infer ontologies \
        of genes downstream of transcription factors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE V3",
        ],
    python_requires='>=3.6',
    install_requires=[
        'matplotlib>=2.0.2',
        'numpy>=1.12.1',
        'scipy>=0.19.1',
        'pandas>=1.0.5'
        ]
)