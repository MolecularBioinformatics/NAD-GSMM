#!/usr/bin/env python
# coding: utf-8

import cofactors

organisms = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
sabio_path = '../generated_data/sabiork_queries/'
brenda_path = '../generated_data/brenda_queries/'
cofactors.brenda_fetch(outpath = brenda_path)
for organism in organisms:
	cofactors.sabio_fetch(organisms = [organism], outpath = sabio_path)
