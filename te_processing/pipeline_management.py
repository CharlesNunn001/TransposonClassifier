#!/usr/bin/env python3

import re
from numpy import *
import pandas as pd
import matplotlib as mp

class PipelineEvaluation:

    def __init__(self, pipe_out, mask_out, out_dir):
        self.pipeline_out = pipe_out
        self.rmasker_out = mask_out
        self.end_dir = out_dir

    def gather_viable_alignments(self):
        ''' Creates a dictionary with masker output high match value lines {entry:{maskername:curatedname}} '''
        self.match_data = []
        with open(self.rmasker_out, 'r') as comparison_entries:
            for match in comparison_entries:
                start, end, diff = match.split('\t')[5:8]
                diff.replace('(','')
                diff.replace(')','')
                diff = int(diff)
                grade = (end - start)/((end + diff) - start)
                if grade >= 0.1 :
                    cur_name = match.split('\t')[8]
                    pipe_name = match.split('\t')[4]
                    lineage = pipe_name.split('#')[0]
                    self.match_data.append([lineage, cur_name])

    def create_curated_list(self):
        ''' Rips curated library into name list '''
        self.curated_names = []
        with open(self.pipeline_out, 'r') as transposon_entries:
            for tag in transposon_entries:
                if tag.startswith('>'):
                    TE_id = tag.split(' ')[0]
                    transposon = TE_id.replace('>', '')
                    self.curated_names.append(transposon)

    def isolate_returned_results(self):
        data_set = []
        with open(self.pipeline_out, 'r') as transposon_entries:
            for tag in transposon_entries:
                if tag.startswith('>'):
                    round_tag = tag.split('\t')[4]
                    trans_id, lineage = round_tag.split('#')
                    if '/' in lineage:
                        order, superfamily = lineage.split('/')
                    else:
                        order = lineage
                        superfamily = 'Unknown'
                    if order == 'LTR':
                        if 'INT' in tag:
                            continue
                    transposon = {'ID': trans_id.replace('>',''), 'Order': order, 'Superfamily': superfamily}
                    data_set.append(transposon)
        transposon_frame = pd.DataFrame(data_set)
        transposon_frame = transposon_frame.set_index('ID')
        transposon_frame = transposon_frame['Order'].replace({'SINE?': 'SINE'}, inplace=True)
        transposon_frame['Counted'] = 0
        transposon_frame['Count'] = 0
        self.pipline_data = transposon_frame

    def find_missing_curations(self):
        temp_names = self.curated_names
        for match in self.match_data:
            if match[1] in temp_names:
                temp_names.remove(match[1])
        return temp_names

    def get_call_counts(self):
        for nameset in self.match_data:
            self.pipline_data[nameset[1]]['Counted'] = 1
            self.pipline_data[nameset[1]]['Count'] += 1
        multiple_matches = self.pipline_data.loc[self.pipline_data['Count'] > 1]
        return multiple_matches

    def make_comparisons(self):
        ''' Manipulate a dataframe from the data_set '''
        fig_seg = self.pipline_data.groupby('Order').count()
        orders = list(fig_seg.index)
        for counter, order in enumerate(orders):
            refined_df = self.pipline_data.loc[self.pipline_data['Order'] == order]
            ds = refined_df.groupby('Counted').count()
            print(ds)
        pass
        #  find called/present for each return and create new frame
        # order | total | called
        # create a plot of these and then save as a file and as a figure
