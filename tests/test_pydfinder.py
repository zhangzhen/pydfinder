"""
Tests for `pydfinder` module.
"""
import mock
from pydfinder import variants
from StringIO import StringIO

class TestVariantFile(object):

    @classmethod
    def setup_class(cls):
        pass

    @mock.patch('pydfinder.variants.open')
    def test_svsim_bedpe_file(self, mock_open):
        expected_d1 = variants.Deletion(
            'DEL0718::1::1',
            variants.GenomeRegion(
                variants.GenomePosition('1', 576974),
                variants.GenomePosition('1', 577575)
            )
        )
        mock_open.return_value = StringIO('1\t576973\t576974\t1\t577574\t577575\tDEL0718::1::1\t255\t+\t+\n')
        with variants.SvsimBedpeFile('any name') as f:
            print f.__class__.__name__
            mock_open.assert_called_once_with('any name', 'r')
            d1 = f.next()
            assert d1 == expected_d1

    @classmethod
    def teardown_class(cls):
        pass
