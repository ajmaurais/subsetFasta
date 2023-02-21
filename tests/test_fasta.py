
import unittest
import re
import os
import tempfile

from .context import subsetFasta, DATA_DIR
from subsetFasta.submodules import fasta as fasta


class TestFasta(unittest.TestCase):

    # path to test fasta file
    TEST_FASTA = f'{DATA_DIR}/test.fasta'

    # uniprot IDs of entries in data/test.fasta
    UNIPROT_ENTRIES = ['P14598', 'P19878', 'Q15080', 'Q01813', 'P17858', 'P11142',
                       'Q9NZV6', 'Q9UJ68', 'P00533', 'Q9H0J9', 'P04406', 'Q96SW2']

    @staticmethod
    def _read_test_file(path):
        fasta_reader = fasta.FastaFile()
        fasta_reader.read(path)
        return fasta_reader

    def test_read_ids(self):
        fasta_reader = self._read_test_file(TestFasta.TEST_FASTA)
        # test individual accessions
        for accession in TestFasta.UNIPROT_ENTRIES:
            self.assertTrue(fasta_reader.id_exists(accession))

        # test iter_ids
        reader_ids = sorted([accession for _, accession in fasta_reader.iter_ids()])
        self.assertListEqual(reader_ids, sorted(TestFasta.UNIPROT_ENTRIES))


    def test_no_dummy_ids(self):
        fasta_reader = self._read_test_file(TestFasta.TEST_FASTA)
        self.assertFalse(fasta_reader.id_exists('DUMMY_ID'))
        with self.assertRaises(KeyError):
            fasta_reader.get_sequence('DUMMY_ID')

    def test_get_sequence(self):
        fasta_reader = self._read_test_file(TestFasta.TEST_FASTA)
        MSRB1_SEQ = 'MSFCSFFGGEVFQNHFEPGVYVCAKCGYELFSSRSKYAHSSPWPAFTETIHADSVAKRPEHNRSEALKVSCGKCGNGLGHEFLNDGPKPGQSRFUIFSSSLKFVPKGKETSASQGH'
        self.assertTrue(fasta_reader.get_sequence('Q9NZV6') == MSRB1_SEQ)

    def test_format_sequence(self):
        fasta_reader = self._read_test_file(TestFasta.TEST_FASTA)
        reader_seqs = {accession: entry[1] for _, accession, entry in fasta_reader.iter_items()}

        # test no max line length
        formated = {accession: fasta.format_sequence(seq) for accession, seq in reader_seqs.items()}
        self.assertDictEqual(formated, reader_seqs)

        # test max line at multiple lengths
        for line_len in (1, 2, 5, 10, 20, 60):
            formated = {accession: fasta.format_sequence(seq, max_line_len=line_len) for accession, seq in reader_seqs.items()}
            for k in formated.keys():
                lines = formated[k].split(sep='\n')
                self.assertTrue(max(len(line) for line in lines) == line_len)
                formated[k] = re.sub(r'\s', '', formated[k])
            self.assertDictEqual(formated, reader_seqs)

    def test_read_formated_sequences(self):
        fasta_reader = self._read_test_file(TestFasta.TEST_FASTA)
        
        for line_len in (None, 1, 5, 10, 60):
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_file_path = os.path.join(temp_dir, 'temp.fasta')
                with open(temp_file_path, 'w') as outF:
                    for i, accession, (description, seq) in fasta_reader.iter_items():
                        seq_f = fasta.format_sequence(seq, line_len)
                        outF.write(f'>sp|{accession}|{description}\n{seq_f}\n')

                test_reader = self._read_test_file(temp_file_path)
                self.assertListEqual(sorted([x for _, x in test_reader.iter_ids()]),
                                     sorted(TestFasta.UNIPROT_ENTRIES))
                for _, accession in test_reader.iter_ids():
                    self.assertEqual(test_reader.get_sequence(accession),
                                     fasta_reader.get_sequence(accession))


if __name__ == '__main__':
    unittest.main(verbosity=2)

