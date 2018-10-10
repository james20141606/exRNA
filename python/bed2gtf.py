#! /usr/bin/env python
import argparse
import sys

def open_file_or_stdout(filename):
    if filename == '-':
        return sys.stdout
    else:
        return open(filename, 'w')

def open_file_or_stdin(filename):
    if filename == '-':
        return sys.stdin
    else:
        return open(filename, 'r')

def bed2gtf(input_file, output_file, source, gene_type):
    bed = open_file_or_stdin(input_file)
    gtf = open_file_or_stdout(output_file)

    lineno = 0
    for line in bed:
        lineno += 1
        c = line.strip().split('\t')
        if len(c) < 12:
            raise ValueError('expect at least 12 columns at line {}'.format(lineno))
        chromStart = int(c[1])
        chromEnd = int(c[2])
        record = [c[0], source, 'gene', chromStart + 1, chromEnd, c[4], c[5], ".",
            'gene_id "{}"; gene_name "{}"; gene_type "{}";"'.format(c[3], c[3], gene_type)]
        gtf.write('\t'.join(map(str, record)) + '\n')
        record = [c[0], source, 'transcript', chromStart + 1, chromEnd, c[4], c[5], ".",
            'gene_id "{}"; transcript_id "{}"; gene_name "{}"; gene_type "{}";"'.format(c[3], c[3], c[3], gene_type)]
        gtf.write('\t'.join(map(str, record)) + '\n')
        blockCount = int(c[9])
        blockSizes = [int(a) for a in c[10].split(',')][:blockCount]
        blockStarts = [int(a) for a in c[11].split(',')][:blockCount]
        exon_number = 0
        for blockSize, blockStart in zip(blockSizes, blockStarts):
            exon_number += 1
            record = [c[0], source, 'exon', chromStart + blockStart + 1, chromStart + blockStart + blockSize, c[4], c[5], ".",
            'gene_id "{}"; transcript_id "{}"; exon_id "{}"; gene_name "{}"; exon_number {}; gene_type "{}";"'.format(c[3],
                 c[3], c[3], c[3], exon_number, gene_type)] 
            gtf.write('\t'.join(map(str, record)) + '\n')
    
    gtf.close()
    bed.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', '-i', type=str, default='-')
    parser.add_argument('--output-file', '-o', type=str, default='-')
    parser.add_argument('--source', type=str, default='BED')
    parser.add_argument('--gene-type', type=str, default='unknown')
    args = parser.parse_args()

    source = args.source
    gene_type = args.gene_type

    try:
        bed2gtf(args.input_file, args.output_file, args.source, args.gene_type)
    except BrokenPipeError:
        pass
    except KeyboardInterrupt:
        pass