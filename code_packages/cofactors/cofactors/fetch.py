#!usr/bin/python
import requests
import os
from zeep import Client
from hashlib import sha256
from pathlib import Path
import pandas as pd
from io import StringIO
from typing import List

def sabio_fetch(
    outpath = './sabiork_queries',
    organisms = []
    ) -> List[str]:
    """
    Queries SabioRK for the specified organisms, writes queries into appropriate files and returns paths to the files.
    :param outpath: Path (folder) into which to save the results.
    :type outpath: str
    :param organisms: List of scientific names of organisms to be queried, e.g. "Homo sapiens"
    :type organisms: List[str]
    :return: List of the files the queries have been written into.
    :rtype: List[str]
    """
    # TODO: only works for first entry in the list, 
    # workaround: use a single element list, e.g. sabio_fetch(['Mus musculus'])
    QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'
    # specify search fields and search terms
    outfiles = []
    outpath = Path(outpath)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    if not organisms:
        organisms = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
    for org in organisms:
        org_split = org.split(' ')
        outfile = outpath / f'sabiork_{org_split[0]}_{org_split[1]}.tsv'
        query_dict = {"Organism": '"{}"'.format(org)}
        # left in for easier modification
        query_string = ' AND '.join(['%s:%s' % (k, v)
                                     for k, v in query_dict.items()])
        # specify output fields and send request
        query = {'fields[]': ['EntryID', 'Organism', 'UniprotID',
                              'ECNumber', 'Parameter'], 'q': query_string}
        request = requests.post(QUERY_URL, params=query)
        request.raise_for_status()
        with open(outfile, 'w') as f:
            f.write(request.text)
        outfiles.append(outfile)

        return outfiles

def brenda_fetch(
    organisms = [], 
    cofactors = [], 
    outpath = './brenda_queries', 
    ) -> List[str]:
    """
    Fetch Km values for given organisms and cofactors from Brenda and save them as tsv files.
    :param organisms: Organisms for which to fetch Km, defaults to []
    :type organisms: List[str], optional
    :param cofactors: Cofactors for which to fetch Km, defaults to []
    :type cofactors: List[str], optional
    :param outpath: Path under which to save TSVs, defaults to './brenda_queries'
    :type outpath: str, optional
    :return: List of paths to all files created
    :rtype: List[str]
    """
    if not cofactors:
        cofactors = ['NAD+', 'NADH', 'NADP+', 'NADPH']    
    if not organisms:
        organisms = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
    email = input('Your Brenda account E-Mail\n')
    password = input('Your Brenda account password\n')
    password = sha256(password.encode("utf-8")).hexdigest()
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir()
    outfiles = []
    for organism in organisms:
        organism_tag = organism.replace(' ', '_')
        for cofactor in cofactors:
            result = _brenda_fetch_single(email, password, organism, cofactor)
            result = _brenda_parse_single(result)
            filename = outpath / f'{cofactor.lower()}_{organism_tag}_brenda.tsv'
            with open(filename, 'w') as f:
                f.write(result)
            outfiles.append(filename)
    return outfiles

        
def _brenda_fetch_single(email:str , password: str, organism: str, cofactor: str):
    """
    Fetch a single set of Kms from Brenda.
    :param email: E-Mail Address of Brenda account
    :type email: str
    :param password: Password of Brenda account (SHA256 HexDigest)
    :type password: str
    :param organism: Organism string (e.g. "Homo sapiens" or "Mus musculus")
    :type organism: str
    :param cofactor: Cofactor string (e.g. "NAD+" or "NADPH")
    :type cofactor: str
    :return: Brenda API's response
    :rtype: Zeep return object (list-like)
    """
    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    client = Client(wsdl)
    parameters = (email, password, "ecNumber*", "kmValue*", "kmValueMaximum*", f"substrate*{cofactor}", 
                  "commentary*", f"organism*{organism}", "ligandStructureId*", "literature*")
    results = client.service.getKmValue(*parameters)
    if not results:
        print(f'List of results is empty for {organism} and {cofactor}. Check correctness of login data.')
    return results

def _brenda_parse_single(result) -> str:
    """
    Parse the response to getKmValue from the Brenda APIs.
    :param result: Single result from Brenda API
    :type result: List-like Zeep return object
    :return: Result parsed into a TSV string
    :rtype: str
    """
    parsed_results = []
    for entry in result:
        km = entry['kmValue']
        literature = entry['literature']
        literature = [str(x) for x in literature]
        literature = ','.join(literature)
        organism = entry['organism']
        ec = entry['ecNumber'] 
        substrate = entry['substrate']
        ls_id = entry['ligandStructureId']
        commentary = entry['commentary']
        if not commentary:
            commentary = ''
        if isinstance(ls_id, int):
            ls_id = str(ls_id)
        entries = [ec, km, substrate, commentary, organism, ls_id, literature]
        entry_str = '\t'.join(entries)
        parsed_results.append(entry_str)
    return '\n'.join(parsed_results)