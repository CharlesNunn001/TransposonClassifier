#!/usr/bin/env python3

import os
import re
from numpy import *
import pandas as pd
import matplotlib as mp

class RGraphGen:

    def __init__(self, inputfile, dirpath):
        self.studies = inputfile
        self.main_directory = dirpath
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

    def process(self):
        class_list = {"DNA": "ClassII", "RNA": "ClassI", "LINE": "ClassI", "LTR": "ClassI", "RC": "ClassII",
                      "rRNA": "ncRNA",
                      "Satellite": "Other", "Simple_repeat": "Other", "SINE": "ClassI", "snRNA": "ncRNA",
                      "tRNA": "ncRNA", "Unknown": "Other", 'Other': "Other", 'ARTEFACT': "Other"}
        for name, study in self.file_pairs.items():
            if study['location'] == None or study['families'] == None:
                continue
            name = name.split(' ')[0]
            df = pd.read_csv(f"{study['location']}/transposon_comparative.tsv", sep='\t')
            df['Order'].replace({'SINE?': 'SINE'}, inplace=True)
            classification = []
            for index, row in df.iterrows():
                classification.append(class_list[row['Order']])
            df['Class'] = classification
            self.construct_main_graphs(df, name, study)
            self.construct_subgraph_pages(df, study, name)

    def construct_subgraph_pages(self, df, study, name):
        colour = ["#ba1a1a", "#404040", "#193791", "#e64040", "#e6e6e6", "#6185f2", "#b09f05", "#5c5c5c", "#ccbb21", "#e8d956",
                  "#f5e873", "#385dc9", "#c4c4c4", "#2d2d2e", "#ba1a1a", "#404040", "#193791", "#e64040", "#e6e6e6", "#6185f2",
                  "#b09f05", "#5c5c5c", "#ccbb21", "#e8d956", "#f5e873", "#385dc9", "#c4c4c4", "#2d2d2e"]
        fig_seg = df.groupby('Order').count()
        orders = list(fig_seg.index)
        cluster = mp.pyplot.figure()
        for counter, order in enumerate(orders):
            refined_df = df.loc[df['Order'] == order]
            ds = refined_df.groupby('Superfamily').count()
            if len(ds.index) == 1 and ds.index[0] == 'Unknown':
                continue
            else:
                ax = cluster.add_subplot(len(orders), 2, (1 + (2 * counter)))
                plot = ds.plot.bar(y='Class', figsize=(25, 25), ax=ax, color=[colour[i] for i in range(len(orders)+1)])
                ax.get_legend().remove()
                ax = cluster.add_subplot(len(orders), 2, (2 + (2 * counter)))
                plot = ds.plot.pie(y='Class', figsize=(30, 30), autopct='', startangle=90, ax=ax,
                                   labels=['' for i in range(len(colour))], colors=colour)
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
        plot = ds.plot.pie(y='Class', figsize=(10, 10), autopct='', labels=['' for i in range(len(indexing)+1)],
                           startangle=90, colors=colours)
        mp.pyplot.xlabel('off')
        mp.pyplot.axis('off')
        mp.pyplot.savefig(f"{study['location']}/pie.png")
        mp.pyplot.clf()
        mp.pyplot.close('all')
