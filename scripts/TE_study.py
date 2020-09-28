#!/usr/bin/env python3

import argparse
from sys import argv

from te_processing.argument_parser import ArgumentParser
from te_processing.rmodeller_strip import dataframeManipulation
from te_processing.family_processing import RGraphGen


def develop(arguments: argparse.Namespace):
    print('Input:', arguments.input)
    print('Output:', arguments.output)
    if arguments.gzip:
        print('Tar Archive:', arguments.gzip)
    else:
        print('Directory:', arguments.directory)
    print('Running TETSV setup...')
    builder = dataframeManipulation(arguments.input, arguments.output, format)
    if arguments.gzip:
        builder.file_process()
    builder.check_outdir_exists()
    builder.dir_process()
    for rmout_file, rmout_loc in builder.file_pairs.items():
        builder.sort_location(rmout_loc['new_loc'])
        builder.isolate_id_history(rmout_loc['fa_loc'], rmout_file)
        builder.save_to_csv(rmout_loc['new_loc'])
    if arguments.gzip:
        builder.clean_up()

def basics(arguments: argparse.Namespace):
    print('Directory:', arguments.dirpath)
    print('Input', arguments.input)
    print('Running family processing...')
    processor = RGraphGen(arguments.input, arguments.dirpath)
    processor.navigate_file_structure()
    processor.process()


def study(arguments: argparse.Namespace):
    print('TEType', arguments.tetype)
    print('Category', arguments.category)
    if arguments.category == 'clade':
        print('Fetching clade list...')
    print('Name', arguments.name)
    print('Input', arguments.input)
    print('Output', arguments.output)


parser = ArgumentParser(develop, basics, study)
args = parser.parse(argv[1:])
if args is not None and args.execute is not None:
    args.execute(args)
