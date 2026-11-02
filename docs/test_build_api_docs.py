"""Deployment checks use rendered files rather than trusting generated URLs."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest

from build_api_docs import build, check_links


class ApiDocumentationTests(unittest.TestCase):
    def test_missing_submodule_target_has_actionable_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'src'
            source.mkdir()
            (source / 'gsw_gibbs.f90').symlink_to('../pkg/GSW-Fortran/gsw_gibbs.f90')
            with self.assertRaisesRegex(FileNotFoundError, 'git submodule update --init --recursive'):
                build(source, root / 'site')

    def test_cross_module_build_and_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'src'
            source.mkdir()
            (source / 'library.F90').write_text(
                'module library\ncontains\nsubroutine helper()\n'
                'end subroutine helper\nend module library\n')
            (source / 'client.F90').write_text(
                'module client\nuse library, only: work => helper\ncontains\n'
                'subroutine run()\ncall work()\ncall unavailable()\n'
                'end subroutine run\nend module client\n')
            output = root / 'site'
            with contextlib.redirect_stdout(io.StringIO()):
                summary = build(source, output)
            self.assertEqual(summary['source_files'], 2)
            self.assertEqual(summary['module_pages'], 2)
            self.assertEqual(summary['html_pages'], 3)
            self.assertEqual(summary['unresolved_references'], 1)
            self.assertGreater(summary['cross_file_links'], 0)
            self.assertIn('library.html#', (output / 'html/client.html').read_text())
            self.assertIn('unavailable', (output / 'reports/unresolved.txt').read_text())
            self.assertEqual((output / 'reports/render.log').read_text(), '')
            self.assertTrue((output / 'rst/index.rst').exists())

    def test_broken_page_and_fragment_links_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'index.html').write_text('<a href="missing.html">Missing page</a>')
            with self.assertRaisesRegex(ValueError, 'missing.html'):
                check_links(root)
            (root / 'missing.html').write_text('<p id="present">Text</p>')
            (root / 'index.html').write_text('<a href="missing.html#absent">Missing anchor</a>')
            with self.assertRaisesRegex(ValueError, 'missing.html#absent'):
                check_links(root)


if __name__ == '__main__':
    unittest.main()
