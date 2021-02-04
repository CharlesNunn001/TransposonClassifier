#!/usr/bin/env python3

import re
from numpy import *
import pandas as pd
import matplotlib as mp

class PipelineEvaluation:

    def __init__(self, pipe_out, mask_out, out_dir, cur_in):
        self.pipeline_out = pipe_out
        self.rmasker_out = mask_out
        self.curated = cur_in
        self.end_dir = out_dir

    @staticmethod
    def convert_to_int(var):
        try:
            var = int(var)
        except:
            var = var.replace('(', '')
            var = var.replace(')','')
            var = int(var)
        return var

    def gather_viable_alignments(self):
        ''' Creates a dictionary with masker output high match value lines {entry:{maskername:curatedname}} '''
        self.match_data = {}
        header_lines = 3
        counter = 0
        with open(self.rmasker_out, 'r') as comparison_entries:
            for match in comparison_entries:
                if header_lines > 0:
                    header_lines -= 1
                    continue
                matchlist = match.split()[1:2] + match.split()[7:10] + match.split()[11:14]
                differ, pipe_name, direct, cur_name, Lstart, Lend, Lremain = matchlist
                lineage = pipe_name.split('#')[0]
                ref_id = lineage + '_' + str(counter)
                counter += 1
                Lstart = self.convert_to_int(Lstart)
                Lend = self.convert_to_int(Lend)
                Lremain = self.convert_to_int(Lremain)
                if cur_name in self.match_data.keys():
                    if direct == '+':
                        self.match_data[cur_name][1][ref_id] = [Lstart, Lend]
                    else:
                        self.match_data[cur_name][1][ref_id] = [Lremain, Lend]
                else:
                    if direct == '+':
                        self.match_data[cur_name] = ((Lend + Lremain), {ref_id: [Lstart, Lend]})
                    else:
                        self.match_data[cur_name] = ((Lend + Lstart), {ref_id: [Lremain, Lend]})

    def create_curated_list(self):
        ''' Rips curated library into name list '''
        self.curated_names = []
        with open(self.curated, 'r') as transposon_entries:
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
                    round_tag = tag.split()[0]
                    trans_id, lineage = round_tag.split('#')
                    if '/' in lineage:
                        order, superfamily = lineage.split('/')
                    else:
                        order = lineage
                        superfamily = 'Unknown'
                    transposon = {'ID': trans_id.replace('>',''), 'Order': order, 'Superfamily': superfamily}
                    data_set.append(transposon)
        transposon_frame = pd.DataFrame(data_set)
        transposon_frame = transposon_frame.set_index('ID')
        transposon_frame['Order'].replace({'SINE?': 'SINE'}, inplace=True)
        transposon_frame['Counted'] = [0 for x in transposon_frame.index]
        transposon_frame['Count'] = [0 for x in transposon_frame.index]
        self.pipline_data = transposon_frame

    @staticmethod
    def find_dist(lower, upper, distance, strip):
        ''' lo/up index False indicates in null (removed) area,
         order: missing, null, partial lower, partial upper, spanning'''
        lo, up = PipelineEvaluation.configure_start_points(lower, upper)
        if lo and up and distance == 1:
            removed = strip[upper] - strip[lower]
        elif not lo and not up and distance == 1:
            return 0, strip
        elif lo and not up and distance > 1:
            removed = strip[lower+1] - strip[lower]
            strip = strip[:lower+1] + strip[upper+1:]
        elif not lo and up and distance > 1:
            removed = strip[upper] - strip[upper-1]
            strip = strip[:lower] + strip[upper:]
        elif not lo and not up and distance > 1:
            removed = strip[upper-1] - strip[lower+1] + 1
            strip = strip[:lower] + strip[upper+1:]
        else:
            return 0, strip
        removed = PipelineEvaluation.clean_edges(removed, strip)
        return removed, strip

    @staticmethod
    def clean_edges(removed, strip):
        '''Removes double edges and adds the remove count for each one'''
        marker = 0
        while marker < len(strip) - 1:
            while strip[marker] == strip[marker + 1]:
                strip.pop(marker + 1)
                strip.pop(marker)
                if marker > 1:
                    marker -= 2
                edge_rem = True
                if len(strip) == 0:
                    break
            marker += 1
        try:
            if edge_rem:
                removed += 1
        except:
            pass
        return removed

    @staticmethod
    def configure_start_points(lower, upper):
        ''' Uses position divisibility to figure out whether the points start in an unmatched or removed area '''
        if (lower % 2) == 0:
            lo = False
        else:
            lo = True
        if upper % 2 == 0:
            up = True
        else:
            up = False
        return lo, up

    def coverage_factors(self): # add unused matches separately
        ''' Restructures the read data from RepeatMasker and finds meta statistics'''
        self.constructed_data = {}
        for ref_id, values in self.match_data.items():
            required_matches = {}
            match_order =[]
            unused_seq = []
            open_range = [1, values[0]]
            best_box = open_range
            while len(values[1]) > 0:
                best_box, largest_removed, master_match, unused_temp = self.calculate_alignment_box(best_box,
                                                                                                    open_range, values)
                unused_seq = unused_seq + unused_temp
                for match in unused_temp:
                    values[1].pop(match)
                try:
                    values[1].pop(master_match)
                except:
                    pass
                open_range = best_box
                match_order.append(master_match)
                if master_match in required_matches.keys():
                    required_matches[master_match][0] = required_matches[master_match][0] + largest_removed
                    required_matches[master_match][1] += 1
                elif master_match != '' :
                    required_matches[master_match] = [largest_removed, 1]
            inverse_coverage = values[0]
            iterbox = iter(best_box)
            for pair in iterbox:
                inverse_coverage = inverse_coverage + (pair - next(iterbox))
            coverage = inverse_coverage/values[0] * 100
            self.constructed_data[ref_id] = (coverage, required_matches, match_order, unused_seq)

                # add matchname to set (coverage, dict: {contig: [total_removed, number_used]}, order used, unused)
                # e.g mat0001: ( 65.1312331,
                #                {jig001: [15, 2], doe004: [300, 6]},
                #                [doe004, doe004, doe004, jig001, doe004, jig001, doe004, doe004],
                #                [doe004, fork006, jig001] )

    def calculate_alignment_box(self, best_box, open_range, values):
        ''' Finds indices of alignments and uses found coverage distances to pick best and null matches for compilation'''
        unused_temp = []
        removed = 0
        largest_removed = 0
        master_match = ''
        for match, locset in values[1].items():
            math_box = open_range + locset
            math_box.sort()
            upper_idx = math_box.index(locset[1])
            lower_idx = math_box[::-1].index(locset[0]) - 1
            dist = upper_idx - lower_idx
            removed, math_box = self.find_dist(lower_idx, upper_idx, dist, math_box)
            if removed == 0:
                unused_temp.append(match)
            if removed > largest_removed:
                best_box = math_box
                largest_removed = removed
                master_match = match
        return best_box, largest_removed, master_match, unused_temp

    def find_missing_curations(self):
        temp_names = self.curated_names
        for match in self.match_data.keys():
            if match in temp_names:
                temp_names.remove(match)
        return temp_names

    def get_call_counts(self):
        for nameset in self.match_data:
            self.pipline_data['Counted'][nameset[0]] = 1
            self.pipline_data['Count'][nameset[0]] += 1
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
            print(order, ds['Order'][0], ds['Order'][1])


    def extra_data(self, missing_tes, multiple_entries):
        multiple_entries.to_csv(f'{self.end_dir}/multiple_entries.tsv', sep='\t')
        with open(f'{self.end_dir}/missed_TEs.txt', 'w') as missedTEs:
            for te in missing_tes:
                missedTEs.write((te + '\n'))

    def formatted_printing(self):
        with open(f'{self.end_dir}/TE_alignment_analysis.txt', 'w') as TEout:
            data1_formatting = []
            for name, data in self.constructed_data.items():
                for matched_te, counters in data[1].items():
                  data1_formatting.append(f'{matched_te}: removed {counters[0]} base pairs across {counters[1]} distinct areas')
                data1 = '\n'.join(data1_formatting)
                data1_formatting = []
                infout = f'''{name}:
Coverage: {data[0]}
Used:
{data1}
Order: {data[2]}
Unused: {data[3]}

'''

                TEout.write(infout)

    def summary_printing(self):
        type_data = {}
        with open('../cele_classification.tsv','r') as class_data:
            for entry in class_data:
                name, accession, TEtype, Subtype = entry.split('\t')
                type_data[accession] = TEtype
        coverage_summary = {}
        for name, data in self.constructed_data.items():
            accession = name.split('.')[0]
            if accession in type_data.keys():
                if type_data[accession] in coverage_summary.keys():
                    coverage_summary[type_data[accession]] = coverage_summary[type_data[accession]] + [data[0]]
                else:
                    coverage_summary[type_data[accession]] = [data[0]]
        with open(f'{self.end_dir}/analysis_summary.txt','w') as cov_out:
            for te_type, coverage_values in coverage_summary.items():
                average_cov = sum(coverage_values)/len(coverage_values)
                counted = len(coverage_values)
                lowest = len([cov for cov in coverage_values if cov <= 25])
                low = len([cov for cov in coverage_values if cov > 25 and cov <= 50])
                high = len([cov for cov in coverage_values if cov > 50 and cov <= 75])
                highest = len([cov for cov in coverage_values if cov > 75])
                write_text = f'''|{te_type} transposable elements:
Counted: {counted} Average: {average_cov}%
Number per 25th percentile margin: {lowest} {low} {high} {highest}

'''
                cov_out.write(write_text)
