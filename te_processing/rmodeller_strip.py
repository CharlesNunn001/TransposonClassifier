#!/usr/bin/env python3

import tarfile
import pandas as pd
import os
import shutil
import re

class dataframeManipulation:

    def __init__(self, input, output, format):
        self.inputfile = input
        self.outputfile = output
        self.file_pairs = {}
        self.data_set = []

    def dir_process(self):
        ''' Find locations for conversion to csv '''
        location = None
        for root, dirs, files in os.walk(self.inputfile):
            for file in files:
                if re.search(f".*-families.fa", file):
                    family_fa = root
                    location = '/'.join(file.split('_')[0:3])
                    self.file_pairs[file] = {'new_loc': location, 'fa_loc': family_fa}


    def sort_location(self, location):
        genus, species, study = location.split('/')
        genus = self.outputfile + '/' + genus
        if not os.path.isdir(genus):
            os.mkdir(genus)
        if not os.path.isdir(genus + '/' + species):
            os.mkdir(genus + '/' + species)
        if not os.path.isdir(genus + '/' + species + '/' + study):
            os.mkdir(genus + '/' + species + '/' + study)

    def file_process(self):
        ''' Build locations for conversion to csv '''
        tar = tarfile.open(self.inputfile, "r:gz")
        tar.extractall()
        tar.close()
        self.inputfile = self.inputfile.replace('.tar.gz', '')

    def clean_up(self):
        ''' Cleans up opened tar archive '''
        shutil.rmtree(self.inputfile)

    def isolate_id_history(self, fa_loc, fa_id):
        self.data_set = []
        fa_file = fa_loc + '/' + fa_id
        with open(fa_file, 'r') as transposon_entries:
            for tag in transposon_entries:
                if tag.startswith('>'):
                    round_tag = tag.split(' ')[0]
                    trans_id, lineage = round_tag.split('#')
                    if '/' in lineage:
                        order, superfamily = lineage.split('/')
                    else:
                        order = lineage
                        superfamily = 'Unknown'
                    transposon = {'ID': trans_id.replace('>',''), 'Order': order, 'Superfamily': superfamily}
                    self.data_set.append(transposon)

    def save_to_csv(self, location):
        transposon_frame = pd.DataFrame(self.data_set)
        transposon_frame = transposon_frame.set_index('ID')
        transposon_frame.to_csv(f'{self.outputfile}/{location}/transposon_comparative.tsv', sep='\t')

    def check_outdir_exists(self):
        if not os.path.isdir(self.outputfile):
            os.mkdir(self.outputfile)
            print(f"""
            WARNING: Output directory not found. A new directory for your TSV has been created at '{self.outputfile}'
""")
