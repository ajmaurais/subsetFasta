
import sys
import os
import argparse

from .submodules.fasta import FastaFile, format_sequence

def readIdsFromFile(fname):
    with open(fname, 'r') as inF:
        lines = [x.strip() for x in inF.readlines()]
    return lines


def writeEntries(fname, entries, multi_line=False):
    with open(fname, 'w') as outF:
        for accession, db, description, seq in entries:
            outF.write('>{}|{}|{}\n{}\n'.format(db, accession, description,
                                                format_sequence(seq, 60) if multi_line else seq))


def main():
    parser = argparse.ArgumentParser('Subset fasta based off a list of uniprot IDs.')
    parser.add_argument('-f', '--file', default=None, help='Obtain ids from a file')
    parser.add_argument('--multiLine', default=False, action='store_true',
                        help='Add line break every 60 characters in protein sequences.')
    parser.add_argument('-o', '--ofname', default=None,
                        help='Output fasta name. Default is <fasta>_subset.fasta')
    parser.add_argument('fasta', help='The fasta file to filter')
    parser.add_argument('ids', nargs='*', help='The IDs to select.')
    args = parser.parse_args()
    
    # make ofname
    ofname = args.ofname
    if ofname is None:
        ofname = '{}_subset.fasta'.format(os.path.splitext(os.path.basename(args.fasta))[0])
    print(ofname)

    # init select_ids
    select_ids = list()
    if args.file is not None:
        select_ids = readIdsFromFile(args.file)
    select_ids += args.ids

    # Only keep the first occurance of duplicate ids
    seen_ids = set()
    unique_select_ids = list()
    for protein_id in select_ids:
        if protein_id not in seen_ids:
            unique_select_ids.append(protein_id)
        seen_ids.add(protein_id)

    fasta = FastaFile()
    fasta.read(args.fasta)

    entries = [(protein_id, *fasta.get_entry(protein_id)) for protein_id in select_ids]
    writeEntries(ofname, entries, multi_line=args.multiLine)

if __name__ == '__main__':
    main()

