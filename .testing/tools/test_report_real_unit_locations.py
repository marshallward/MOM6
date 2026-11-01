import unittest
from pathlib import Path
import tempfile

import flint

from report_real_units import real_unit_issues


class TestRealUnitLocations(unittest.TestCase):
    def issues_for(self, declaration):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.F90'
            path.write_text('subroutine demo()\n' + declaration + '\nend subroutine demo\n')
            return real_unit_issues(flint.parse(str(path)).sources[0])

    def test_units_do_not_leak_between_continued_variables(self):
        issues = self.issues_for('  real :: a, & !< Length [m]\n    b !< Undocumented')
        self.assertEqual([issue.name for issue in issues], ['b'])
        self.assertEqual(issues[0].comments, ('!< Undocumented',))

    def test_competing_doxygen_only_for_entity_with_backward_comment(self):
        issues = self.issues_for('!> Shared [m]\n  real :: a, &\n    b !< Specific [m]')
        self.assertEqual([(issue.name, issue.message) for issue in issues],
                         [('b', 'competing Doxygen comments')])

    def test_shared_and_following_units(self):
        for declaration in ('  real :: a, b ! Length [m]',
                            '  real :: a\n  ! Length [m]'):
            self.assertEqual(self.issues_for(declaration), [])

    def test_shared_units_on_continuation_line(self):
        issues = self.issues_for(
            '  real :: &\n'
            '    dvdx, dudy, & ! Circulation contributions [L2 T-1 ~> m2 s-1]\n'
            '    dvSdx, duSdy, & ! Stokes drift contributions [L2 T-1 ~> m2 s-1]\n'
            '    undocumented ! No units'
        )
        self.assertEqual([issue.name for issue in issues], ['undocumented'])

    def test_double_precision_variables_are_real(self):
        issues = self.issues_for('  double precision a, b')
        self.assertEqual([issue.name for issue in issues], ['a', 'b'])

    def test_unrelated_following_documentation_is_not_units(self):
        issues = self.issues_for('  real :: a\n\n  ! Unrelated [m]')
        self.assertEqual([issue.name for issue in issues], ['a'])

    def test_declarator_location_not_earlier_same_spelling(self):
        issues = self.issues_for('real, dimension(x) :: &\n  x, y')
        self.assertEqual([(issue.name, issue.line_number) for issue in issues], [('x', 3), ('y', 3)])

    def test_standard_declarations_use_entity_locations(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.F90'
            path.write_text('subroutine demo()\n'
                            '  real, dimension(x) :: &\n'
                            '    x, y\n'
                            '  real a, b(2)\n'
                            '  real :: documented !< Length [m]\n'
                            'end subroutine demo\n')
            source = flint.parse(str(path)).sources[0]
            issues = real_unit_issues(source)
            self.assertEqual([(issue.name, issue.line_number) for issue in issues],
                             [('x', 3), ('y', 3), ('a', 4), ('b', 4)])


if __name__ == '__main__':
    unittest.main()
