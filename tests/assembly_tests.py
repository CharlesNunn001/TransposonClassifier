import os
import unittest

from te_processing.pipeline_management import PipelineEvaluation
from unittest.mock import patch, MagicMock


class Test(unittest.TestCase):

    def setUp(self):
        tester = PipelineEvaluation('../te_processing/caenorhabditis_elegans_prjna13758_core_15_98_276-families.fa',
                                    '../te_processing/caenorhabditis_elegans_prjna13758_core_15_98_276-families.fa.out',
                                    '../te_processing/', '../te_processing/families.fa')
        self.pipeline_tester = tester

    #Length calc reference '-' = unmatched
    # |   ----------------                    -----     -----------|  Active layer
    # |     -------  --------    ------     ---------        ------|  Adding Layer
    #         1          2          3           4              5
    # 1. Missing piece    2. Partial Piece    3. Redundant piece   4. Spanning Piece   5. Edge Piece

    def test_lengthcalc_missing_piece(self):
        expected = 4, [1, 6, 10, 20]
        distance = 1
        lower = 1
        upper = 2
        strip = [1, 6, 10, 20]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_partial_lo_piece(self):
        expected = 2, [1, 4, 12, 20]
        distance = 2
        lower = 1
        upper = 3
        strip = [1, 4, 6, 10, 12, 20]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_partial_up_piece(self):
        expected = 2, [1, 6, 14, 20]
        distance = 2
        lower = 2
        upper = 4
        strip = [1, 6, 10, 12, 14, 20]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_redundant_piece(self):
        expected = 0, [1, 4, 6, 10, 12, 20]
        distance = 1
        lower = 2
        upper = 3
        strip = [1, 4, 6, 10, 12, 20]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_spanning_piece(self):
        expected = 5, [1, 4, 12, 20]
        distance = 3
        lower = 2
        upper = 5
        strip = [1, 4, 5, 6, 10, 11, 12, 20]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_edge_piece(self):
        expected = 3, [5,10]
        distance = 1
        lower = 1
        upper = 2
        strip = [1,1,3,3,5,10]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_lengthcalc_innerjoin_piece(self):
        expected = 3, [1, 2, 7, 10]
        distance = 1
        lower = 3
        upper = 4
        strip = [1, 2, 4, 4, 6, 6, 7, 10]
        actual = self.pipeline_tester.find_dist(lower,upper,distance,strip)
        self.assertEqual(expected, actual)

    def test_parentheses_removed(self):
        parenthesis = '(123)'
        expected = 123
        actual = self.pipeline_tester.convert_to_int(parenthesis)
        self.assertEqual(expected, actual)

    def test_no_fail_when_no_parentheses_removed(self):
        parenthesis = '123'
        expected = 123
        actual = self.pipeline_tester.convert_to_int(parenthesis)
        self.assertEqual(expected, actual)

    def test_clean_double_edges(self):
        edged_piece = [1,1,5,10]
        removed = 5
        expected = 6
        actual = self.pipeline_tester.clean_edges(removed, edged_piece)
        self.assertEqual(expected, actual)

    def test_start_point_config_low_in_up_out(self):
        lower_position = 1
        upper_position = 1
        expected = True, False
        actual = self.pipeline_tester.configure_start_points(lower_position, upper_position)
        self.assertEqual(expected, actual)

    def test_start_point_config_low_out_up_in(self):
        lower_position = 2
        upper_position = 2
        expected = False, True
        actual = self.pipeline_tester.configure_start_points(lower_position, upper_position)
        self.assertEqual(expected, actual)

    def test_new_best_box(self):
        lower_position = 2
        upper_position = 2
        expected = False, True
        actual = self.pipeline_tester.configure_start_points(lower_position, upper_position)
        self.assertEqual(expected, actual)