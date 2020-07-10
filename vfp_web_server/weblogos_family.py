# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:10:05 2020

@author: pedro
"""

from weblogo import *
from Bio.Align.Applications import ClustalOmegaCommandline

def WebLogo(family):
    in_file = family + ".fasta"
    out_file = family + ".clustal_num"
    """
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=False)
    print(clustalomega_cline)
    os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' + str(clustalomega_cline) + ' --force')
    
    file_out = open(out_file, "r")
    # seqs_aligned = file_out.readlines()
    # return_data = {'data': seqs_aligned}
    """

    file_out = open(out_file, "r")
    seqs = read_seq_data(file_out)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "VFP WEBSERVER"
    logoformat = LogoFormat(logodata, logooptions)
    eps = txt_formatter(logodata, logoformat)

    weblogo_file = family + "_weblogo.txt"
    weblogo = open(weblogo_file, "w")
    data_weblogo = str(eps)[2:len(str(eps))-1].replace('\\n', '\n').replace('\\t', '\t')
    weblogo.write(data_weblogo)
    weblogo.close()

if __name__ == '__main__':
    for i in ['Arenaviridae', 'Coronaviridae', 'Filoviridae' , 'Flaviviridae',
              'Herpesviridae', 'Orthomyxoviridae', 'Paramyxoviridae', 'Peribunyaviridae', 'Phenuiviridae', 'Pneumoviridae',
              'Retroviridae', 'Rhabdoviridae', 'Togaviridae']:
        WebLogo(i)