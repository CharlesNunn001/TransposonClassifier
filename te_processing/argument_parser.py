import argparse


class ArgumentParser:

    def __init__(self, setup=None, basic=None, study=None):
        self.setup = setup
        self.basic = basic
        self.study = study
        self.parser = self._build_parser()
        pass

    def parse(self, args):
        try:
            return self.parser.parse_args(args)
        except SystemExit:
            return None

    def _build_parser(self):
        parser = argparse.ArgumentParser(prog='TransposonClassifier')
        subparsers = parser.add_subparsers(help='sub-command help')
        self._set_up_parser(subparsers.add_parser('setup', help='Sets up TSV for further steps'))
        self._basic_parser(subparsers.add_parser('build', help='Basic graphing on all data in species list'))
        self._study_parser(subparsers.add_parser('study', help='Print the command to import the external data'))
        return parser

    def _study_parser(self, study):
        study.add_argument('-t', '--tetype', required=True, help='The transposable elements to be looked at')
        study.add_argument('-c', '--category', choices=['clade', 'genus', 'species'], required=True,
                           help='The clade, genus or species to be studied')
        study.add_argument('-n', '--name', required=True, help='The name of category')
        study.add_argument('-i', '--input', required=True, help='Base directory for tsv files')
        study.add_argument('-o', '--output', required=True, help='Directory for command file')
        study.set_defaults(execute=self.study)

    def _basic_parser(self, basics):
        basics.add_argument('-d', '--dirpath', required=True, help='The location of the top directory for the studies')
        basics.add_argument('-i', '--input', required=True, help='The input studies and categories file')
        basics.set_defaults(execute=self.basic)

    def _set_up_parser(self, set_up):
        set_up.add_argument('-i', '--input', required=True, help='Input file or top level directory')
        set_up.add_argument('-o', '--output', required=True, help='Top level output directory')
        group = set_up.add_mutually_exclusive_group(required=True)
        group.add_argument('-d', '--directory', help='Receive fastas from a pre-existing directory structure' ,
                           action='store_true')
        group.add_argument('-gz', '--gzip', help='Receive fastas in a surface level tar.gzip', action='store_true')
        set_up.set_defaults(execute=self.setup)