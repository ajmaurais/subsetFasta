
import sys
import os
import argparse

from submodules.fasta import FastaFile

def readIdsFromFile(fname):
    with open(fname, 'r') as inF:
        lines = [x.strip() for x in inF.readlines()]
    return set(lines)


def writeEntries(fname, entries):
    with open(fname, 'w') as outF:
        for entry in entries:
            outF.write('>sp|{}|{}\n{}\n'.format(*entry))


def main():
    parser = argparse.ArgumentParser('Subset fasta based off a list of uniprot IDs.')
    parser.add_argument('-f', '--file', default=None, help='Obtain ids from a file')
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
    writeEntries(ofname, entries)

if __name__ == '__main__':
    main()

