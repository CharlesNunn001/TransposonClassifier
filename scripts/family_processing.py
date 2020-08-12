#!/usr/bin/env python3

import argparse
import os
import re
from numpy import *
import pandas as pd
import matplotlib as mp

class RGraphGen:

    def __init__(self, args):
        self.studies = args.input
        self.main_directory = args.dirpath
        self.file_pairs = {}
        self.types = ['I', 'C', 'III', 'IV', 'V', 'Trematoda', 'Cestoda', 'Monogenea', 'Free-living-flatworm']

    def navigate_file_structure(self):
        with open(self.studies, 'r') as std_file:
            for study in std_file:
                print(study)
                path = '/'.join(((study.split('\t')[0]).split('_'))[0:3])
                spec_type = study.split('\t')[1].replace('\n','')
                RM_file = None
                location = None
                for root, dirs, files in os.walk(self.main_directory + '/' + path):
                    for file in files:
                        if re.search(f".*-families.fa", file):
                            RM_file = root + '/' + file
                            location = root
                            break
                self.file_pairs[study] = {'location': location, 'families': RM_file, 'type': spec_type}

    def r_setup(self):
        for category in self.types:
            if not os.path.isdir(self.main_directory + '/' + category):
                os.mkdir(self.main_directory + '/' + category)
                print(f"New Directory: {category} has been created at '{self.main_directory}'")

    def process(self):
        class_list = {"DNA": "ClassII", "RNA": "ClassI", "LINE": "ClassI", "LTR": "ClassI", "RC": "ClassII",
                      "rRNA": "ncRNA",
                      "Satellite": "Other", "Simple_repeat": "Other", "SINE": "ClassI", "snRNA": "ncRNA",
                      "tRNA": "ncRNA", "Unknown": "Other", 'Other': "Other", 'ARTEFACT': "Other"}
        for name, study in self.file_pairs.items():
            if study['location'] == None or study['families'] == None:
                continue
            name = name.split('\t')[0]
            df = pd.read_csv(f"{study['location']}/transposon_comparative.tsv", sep='\t')
            df['Order'].replace({'SINE?': 'SINE'}, inplace=True)
            classification = []
            for index, row in df.iterrows():
                classification.append(class_list[row['Order']])
            df['Class'] = classification
            self.construct_main_graphs(df, name, study)
            self.construct_subgraph_pages(df, study)

    def construct_subgraph_pages(self, df, study):
        colour = ["#ba1a1a", "#404040", "#193791", "#e64040", "#e6e6e6", "#6185f2", "#b09f05", "#5c5c5c", "#ccbb21", "#e8d956",
                  "#f5e873", "#385dc9", "#c4c4c4", "#2d2d2e"]
        fig_seg = df.groupby('Order').count()
        orders = list(fig_seg.index)
        cluster = mp.pyplot.figure()
        for counter, order in enumerate(orders):
            refined_df = df.loc[df['Order'] == order]
            ds = refined_df.groupby('Superfamily').count()
            ax = cluster.add_subplot(len(orders), 2, (1 + (2 * counter)))
            plot = ds.plot.bar(y='Class', figsize=(10, 10), ax=ax, color=[colour[i] for i in range(len(orders)+1)])
            ax.get_legend().remove()
            ax = cluster.add_subplot(len(orders), 2, (2 + (2 * counter)))
            plot = ds.plot.pie(y='Class', figsize=(15, 15), autopct='', startangle=90, ax=ax,
                               labels=['' for i in range(len(orders)+1)], colors=colour)
            ax.get_legend().remove()
            mp.pyplot.title(f'{order}')
        mp.pyplot.tight_layout()
        mp.pyplot.savefig(f"{study['location']}/detailed_plots.png")
        mp.pyplot.clf()
        mp.pyplot.close('all')


    def construct_main_graphs(self, df, name, study):
        class_colour = {"DNA": "#ba1a1a", "RC": "#e64040", "RNA": "#b09f05", "LINE": "#ccbb21", "LTR": "#e8d956",
                        "SINE": "#f5e873",
                        "rRNA": "#6185f2", "snRNA": "#193791", "tRNA": "#385dc9",
                        "Satellite": "#c4c4c4", "Simple_repeat": "#2d2d2e",
                        "Unknown": "#5c5c5c", 'Other': "#404040", 'ARTEFACT': "#e6e6e6"}
        df = df.sort_values(by=['Class'])
        ds = df.groupby('Order').count()
        indexing = list(ds.index)
        colours = []
        for feature in indexing:
            colours.append(class_colour[feature])
        plot = ds.plot.bar(y='Class', color=colours)
        mp.pyplot.savefig(f"{study['location']}/bar.png")
        mp.pyplot.clf()
        mp.pyplot.savefig(f"{self.main_directory}/{study['type']}/{name}_bar.png")
        mp.pyplot.clf()
        plot = ds.plot.pie(y='Class', figsize=(10, 10), autopct='', labels=['' for i in range(len(indexing)+1)],
                           startangle=90, colors=colours)
        mp.pyplot.xlabel('off')
        mp.pyplot.axis('off')
        mp.pyplot.savefig(f"{study['location']}/pie.png")
        mp.pyplot.clf()
        mp.pyplot.savefig(f"{self.main_directory}/{study['type']}/{name}_pie.png")
        mp.pyplot.clf()
        mp.pyplot.close('all')

    def strip_fasta_files(self):
        for rmout_file in self.file_pairs.values():
            if rmout_file['location'] == None or rmout_file['families'] == None:
                continue
            head_read = dataframeManipulation(rmout_file['families'], rmout_file['location'])
            head_read.isolate_id_history()
            head_read.check_outdir_exists()
            head_read.save_to_csv()

class dataframeManipulation:

    def __init__(self, inputfile, outputfile):
        self.inputfile = inputfile
        self.outputfile = outputfile
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
                    transposon = {'ID': trans_id.replace('>', ''), 'Order': order, 'Superfamily': superfamily}
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
parser.add_argument('-i', '--input', required=True, help='The input studies and categories file')
parser.add_argument('-d', '--dirpath', required=True, help='The location of the top directory for the studies')
args = parser.parse_args()
processor = RGraphGen(args)
processor.navigate_file_structure()
processor.r_setup()
processor.strip_fasta_files()
processor.process()
