
import sys
import os
import argparse

from .submodules.fasta import FastaFile, format_sequence

def readIdsFromFile(fname):
    with open(fname, 'r') as inF:
        lines = [x.strip() for x in inF.readlines()]
    return set(lines)


def writeEntries(fname, entries, multi_line=False):
    with open(fname, 'w') as outF:
        for accession, description, seq in entries:
            outF.write('>sp|{}|{}\n{}\n'.format(accession, description,
                                                format_sequence(seq, 60) if multi_line else seq))


def main():
    parser = argparse.ArgumentParser('Subset fasta based off a list of uniprot IDs.')
    parser.add_argument('-f', '--file', default=None, help='Obtain ids from a file')
    parser.add_argument('--multiLine', default=False, action='store_true',
                        help='Add line break every 60 characters in protein sequences.')
    parser.add_argument('fasta', help='The fasta file to filter')
    parser.add_argument('ids', nargs='*', help='The IDs to select.')
    args = parser.parse_args()
    
    # make ofname
    ofname = '{}_subset.fasta'.format(os.path.splitext(os.path.basename(args.fasta))[0])
    print(ofname)

    # init select_ids
    select_ids = set()
    if args.file is not None:
        select_ids = readIdsFromFile(args.file)
    select_ids.update(set(args.ids))
    fasta = FastaFile()
    fasta.read(args.fasta)

    entries = [(protein, *fasta.get_entry(protein)) for protein in select_ids]
    writeEntries(ofname, entries, multi_line=args.multiLine)

if __name__ == '__main__':
    main()

