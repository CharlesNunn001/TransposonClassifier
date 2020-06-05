#!/usr/bin/env python3

import argparse
import pandas as pd
import os

class dataframeManipulation:

    def __init__(self, args: argparse.Namespace):
        self.inputfile = args.input
        self.outputfile = args.output
        self.data_set = []

    def isolate_id_history(self):
        with open(self.inputfile, 'r') as transposon_entries:
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

    def save_to_csv(self):
        transposon_frame = pd.DataFrame(self.data_set)
        transposon_frame = transposon_frame.set_index('ID')
        transposon_frame.to_csv(f'{self.outputfile}/transposon_comparative.tsv', sep='\t')

    def check_outdir_exists(self):
        if not os.path.isdir(self.outputfile):
            os.mkdir(self.outputfile)
            print(f"""
            WARNING: Output directory not found. A new directory for your TSV has been created at '{self.outputfile}'
""")

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='The input fasta protein file')
parser.add_argument('-o', '--output', required=True, help='The output directory for the tsv file')
args = parser.parse_args()
head_read = dataframeManipulation(args)
head_read.isolate_id_history()
head_read.check_outdir_exists()
head_read.save_to_csv()
