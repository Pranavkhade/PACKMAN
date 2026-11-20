import io
import math
import unittest
import warnings
from pathlib import Path

from ...apps.predict_hinge import predict_hinge
from ...molecule import molecule


PAPER_HINGE_TOLERANCE = 3
PAPER_BENCHMARK_CASES = (
    ('1K20', '1k20.pdb', 'pdb', 'A', 2.6, ((180, 211),)),
    ('1K23', '1k23.pdb', 'pdb', 'A', 3.0, ((185, 195),)),
    ('1PRW', '1prw.pdb', 'pdb', 'A', 2.8, ((71, 96),)),
    ('1EXR', '1exr.cif', 'cif', 'A', 2.8, ((64, 84),)),
    ('3HVP', '3hvp.pdb', 'pdb', 'A', 3.4, ((21, 33), (82, 90))),
    ('1SSP', '1ssp.pdb', 'pdb', 'E', 2.8, ((137, 144), (190, 205), (151, 161))),
    ('1AKZ', '1akz.pdb', 'pdb', 'A', 2.8, ((136, 143), (189, 205), (152, 162))),
    ('2DRI', '2dri.pdb', 'pdb', 'A', 2.8, ((102, 110), (129, 145), (233, 240))),
    ('1URP', '1urp.pdb', 'pdb', 'A', 3.0, ((100, 112), (135, 148), (232, 244))),
)


def _find_hinge_match(expected_range, observed_ranges, tolerance):
    for index, observed_range in enumerate(observed_ranges):
        if all(
            abs(expected - observed) <= tolerance
            for expected, observed in zip(expected_range, observed_range)
        ):
            return (index,)

    ordered_ranges = sorted(enumerate(observed_ranges), key=lambda item: item[1][0])
    merged_ranges = []
    for index, (start, end) in ordered_ranges:
        if merged_ranges and start - merged_ranges[-1][1] - 1 <= tolerance:
            merged_ranges[-1][1] = max(merged_ranges[-1][1], end)
            merged_ranges[-1][2].append(index)
        else:
            merged_ranges.append([start, end, [index]])

    for start, end, indices in merged_ranges:
        if all(
            abs(expected - observed) <= tolerance
            for expected, observed in zip(expected_range, (start, end))
        ):
            return tuple(indices)
    return ()


class TestPredictHinge(unittest.TestCase):

    def test_paper_benchmarks_match_expected_hinges_within_three_residues(self):
        data_dir = Path(__file__).resolve().parents[1] / 'data'

        for pdb_id, filename, filetype, chain_id, alpha, expected_ranges in PAPER_BENCHMARK_CASES:
            with self.subTest(pdb_id=pdb_id):
                structure = molecule.load_structure(str(data_dir / filename), ftype=filetype)
                chain = structure[0][chain_id]
                backbone = [
                    atom
                    for group in chain.get_backbone()
                    for atom in group
                    if atom is not None
                ]
                self.assertTrue(backbone, 'The benchmark chain has no backbone atoms.')

                output = io.StringIO()
                alpha_shape = predict_hinge(
                    backbone,
                    output,
                    Alpha=alpha,
                    filename=filename,
                )
                self.assertTrue(alpha_shape)

                remaining_ranges = [
                    (hinge.get_elements()[0].get_id(), hinge.get_elements()[-1].get_id())
                    for hinge in chain.get_hinges()
                ]
                for expected_range in expected_ranges:
                    matched_indices = _find_hinge_match(
                        expected_range,
                        remaining_ranges,
                        PAPER_HINGE_TOLERANCE,
                    )
                    self.assertTrue(
                        matched_indices,
                        '{} expected hinge {} within {} residues; predicted {}'.format(
                            pdb_id,
                            expected_range,
                            PAPER_HINGE_TOLERANCE,
                            remaining_ranges,
                        ),
                    )
                    for index in sorted(matched_indices, reverse=True):
                        remaining_ranges.pop(index)

    def test_calmodulin_prediction_overlaps_reference_hinge_region(self):
        cases = (
            ('1exr.cif', 'cif', 584),
            ('1prw.pdb', 'pdb', 588),
        )
        reference_hinge = (71, 90)

        for filename, filetype, expected_backbone_count in cases:
            with self.subTest(filename=filename):
                structure = molecule.load_structure(
                    'packman/tests/data/' + filename,
                    ftype=filetype,
                )
                self.assertEqual(len(list(structure)), 1)
                chains = list(structure[0].get_chains())
                self.assertEqual(tuple(chain.get_id() for chain in chains), ('A',))
                chain = chains[0]
                backbone = [
                    atom
                    for group in chain.get_backbone()
                    for atom in group
                    if atom is not None
                ]
                self.assertEqual(len(backbone), expected_backbone_count)

                output = io.StringIO()
                alpha_shape = predict_hinge(
                    backbone,
                    output,
                    Alpha=2.8,
                    filename=filename,
                )
                self.assertTrue(alpha_shape)

                hinges = chain.get_hinges()
                hinge_ranges = tuple(
                    (hinge.get_elements()[0].get_id(), hinge.get_elements()[-1].get_id())
                    for hinge in hinges
                )
                overlaps = tuple(
                    max(0, min(end, reference_hinge[1]) - max(start, reference_hinge[0]) + 1)
                    for start, end in hinge_ranges
                )
                self.assertGreaterEqual(max(overlaps, default=0), 10)

                report_lines = output.getvalue().splitlines()
                self.assertIn(
                    'Filename= ' + filename + '\t| Chain(s)= A\t| AlphaValue= 2.8\t| MinimumHingeLength= 5\t| EccentricityClusters= 4',
                    report_lines,
                )
                self.assertEqual(len(hinges), sum(line.startswith('Hinge #') for line in report_lines))
                for hinge_index, (hinge, hinge_range) in enumerate(zip(hinges, hinge_ranges), start=1):
                    self.assertEqual(hinge.get_id(), hinge_index - 1)
                    self.assertGreater(hinge_range[1] - hinge_range[0] + 1, 5)
                    self.assertIn('Hinge #' + str(hinge_index), report_lines)

    def test_prediction_is_reproducible_and_matches_expected_output(self):
        runs = []
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in range(2):
                structure = molecule.load_structure('packman/tests/data/4hla.cif', ftype='cif')
                chain = structure[0]['A']
                backbone = [
                    atom
                    for group in chain.get_backbone()
                    for atom in group
                    if atom is not None
                ]
                output = io.StringIO()
                alpha_shape = predict_hinge(
                    backbone,
                    output,
                    Alpha=2.8,
                    filename='4hla.cif',
                )
                hinges = chain.get_hinges()
                runs.append({
                    'backbone': backbone,
                    'alpha_shape': alpha_shape,
                    'hinges': hinges,
                    'domains': tuple(
                        (residue.get_id(), residue.get_domain_id())
                        for residue in chain.get_residues()
                    ),
                    'report': output.getvalue(),
                })

        self.assertFalse(any(issubclass(item.category, ResourceWarning) for item in caught))
        self.assertFalse(
            any(
                'KMeans is known to have a memory leak on Windows with MKL' in str(item.message)
                for item in caught
            )
        )
        self.assertFalse(
            any(
                'default behavior of `mode`' in str(item.message)
                for item in caught
            )
        )
        self.assertEqual(len(runs[0]['backbone']), 396)
        self.assertEqual(len(runs[1]['backbone']), 396)
        self.assertTrue(all(len(tetrahedron) == 4 for tetrahedron in runs[0]['alpha_shape']))
        self.assertEqual(len(runs[0]['alpha_shape']), 761)
        input_atom_ids = {atom.get_id() for atom in runs[0]['backbone']}
        self.assertTrue(
            all(
                atom.get_id() in input_atom_ids
                for tetrahedron in runs[0]['alpha_shape']
                for atom in tetrahedron
            )
        )
        alpha_shape_signatures = tuple(
            tuple(sorted(tuple(sorted(atom.get_id() for atom in tetrahedron)) for tetrahedron in run['alpha_shape']))
            for run in runs
        )
        self.assertEqual(alpha_shape_signatures[0], alpha_shape_signatures[1])

        expected_hinge_ranges = (
            tuple(range(21, 35)),
            tuple(range(57, 63)),
            tuple(range(72, 78)),
            tuple(range(82, 92)),
        )
        expected_hinges = tuple(
            (hinge.get_id(), tuple(residue.get_id() for residue in hinge.get_elements()))
            for hinge in runs[0]['hinges']
        )
        self.assertEqual(
            expected_hinges,
            tuple((hinge_id, residue_ids) for hinge_id, residue_ids in enumerate(expected_hinge_ranges)),
        )
        self.assertEqual(
            tuple(
                (hinge.get_id(), tuple(residue.get_id() for residue in hinge.get_elements()))
                for hinge in runs[1]['hinges']
            ),
            expected_hinges,
        )

        for hinge in runs[0]['hinges']:
            self.assertEqual(hinge.get_alpha_value(), 2.8)
            self.assertTrue(math.isfinite(hinge.get_pvalue()))
            self.assertGreaterEqual(hinge.get_pvalue(), 0.0)
            self.assertLessEqual(hinge.get_pvalue(), 1.0)
            stats = hinge.get_stats()
            self.assertEqual(stats[0], ['', 'N', 'Min', 'Max', 'Mean', 'Mode', 'Median', 'STDDev'])
            self.assertEqual([row[0] for row in stats[1:]], ['Total', 'Hinge', 'NonHinge'])
            self.assertEqual([row[1] for row in stats[1:]], [396, 4 * len(hinge.get_elements()), 396 - 4 * len(hinge.get_elements())])

        expected_domains = (
            *((residue_id, 'DM0') for residue_id in range(1, 21)),
            *((residue_id, 'FL0') for residue_id in range(21, 35)),
            *((residue_id, 'DM1') for residue_id in range(35, 57)),
            *((residue_id, 'FL1') for residue_id in range(57, 63)),
            *((residue_id, 'DM2') for residue_id in range(63, 72)),
            *((residue_id, 'FL2') for residue_id in range(72, 78)),
            *((residue_id, 'DM3') for residue_id in range(78, 82)),
            *((residue_id, 'FL3') for residue_id in range(82, 92)),
            *((residue_id, 'DM4') for residue_id in range(92, 100)),
        )
        self.assertEqual(runs[0]['domains'], expected_domains)
        self.assertEqual(runs[1]['domains'], expected_domains)

        expected_report_lines = (
            'Filename= 4hla.cif\t| Chain(s)= A\t| AlphaValue= 2.8\t| MinimumHingeLength= 5\t| EccentricityClusters= 4',
            'Hinge #1',
            'Residues: GLU-21 to GLU-34',
            'Hinge #2',
            'Residues: ARG-57 to ILE-62',
            'Hinge #3',
            'Residues: ILE-72 to VAL-77',
            'Hinge #4',
            'Residues: VAL-82 to THR-91',
        )
        report_lines = runs[0]['report'].splitlines()
        for expected_line in expected_report_lines:
            self.assertIn(expected_line, report_lines)
        self.assertEqual(runs[0]['report'], runs[1]['report'])


if __name__ == '__main__':
    unittest.main()
