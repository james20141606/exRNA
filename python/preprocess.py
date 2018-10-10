#! /usr/bin/env python
from __future__ import print_function
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def transcript_counts(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout
    from collections import defaultdict

    logger.info('read input transcript BAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    counts = defaultdict(int)
    for read in sam.fetch():
        counts[read.reference_name] += 1
    
    logger.info('create output file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        for key, val in counts.items():
            fout.write('{}\t{}\n'.format(key, val))

@command_handler
def gtf_to_transcript_table(args):
    from ioutils import open_file_or_stdin, open_file_or_stdout
    from collections import OrderedDict

    feature = args.feature
    default_transcript_type = args.transcript_type
    default_gene_type = args.gene_type

    fout = open_file_or_stdout(args.output_file)
    with open_file_or_stdin(args.input_file) as fin:
        transcripts = OrderedDict()
        for line in fin:
            c = line.strip().split('\t')
            if c[0].startswith('#'):
                continue
            if c[2] != feature:
                continue
            attrs = {}
            for a in c[8].split(';')[:-1]:
                a = a.strip()
                i = a.find(' ')
                key = a[:i]
                val = a[(i + 1):].strip('"')
                attrs[key] = val
            if 'transcript_name' not in attrs:
                attrs['transcript_name'] = attrs['transcript_id']
            if 'gene_name' not in attrs:
                attrs['gene_name'] = attrs['gene_id']
            if default_transcript_type is not None:
                attrs['transcript_type'] = default_transcript_type
            else:
                if 'transcript_type' not in attrs:
                    attrs['transcript_type'] = 'unknown'
            if default_gene_type is not None:
                attrs['gene_type'] = default_gene_type
            else:
                if 'gene_type' not in attrs:
                    attrs['gene_type'] = 'unknown'
            exon = [c[0], int(c[3]) - 1, int(c[4]), attrs['gene_id'], 0, c[6],
                attrs['gene_id'], attrs['transcript_id'], 
                attrs['gene_name'], attrs['transcript_name'],
                attrs['gene_type'], attrs['transcript_type']]
            transcript = transcripts.get(attrs['transcript_id'])
            if transcript is None:
                transcripts[attrs['transcript_id']] = exon
            else:
                if c[2] == 'exon':
                    transcript[1] = min(transcript[1], exon[1])
                    transcript[2] = max(transcript[2], exon[2])
        header = ['chrom', 'start', 'end', 'name', 'score', 'strand',
            'gene_id', 'transcript_id', 
            'gene_name', 'transcript_name',
            'gene_type', 'transcript_type'
        ]
        print('\t'.join(header), file=fout)
        for transcript in transcripts.values():
            print('\t'.join(str(a) for a in transcript), file=fout)
    fout.close()

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Preprocessing module')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('transcript_counts', 
        help='count reads that belongs to a transcript')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='input transcript BAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output transcript counts file')
    
    parser = subparsers.add_parser('gtf_to_transcript_table')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input GTF file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output table file')
    parser.add_argument('--feature', type=str,
        help='feature to use in input GTF file (Column 3)')
    parser.add_argument('--gene-type', type=str,
        help='gene type to set if "gene_type" attribute is not available in GTF file')
    parser.add_argument('--transcript-type', type=str, default='unknown',
        help='gene type to set if "transcript_type" attribute is not available in GTF file')

    args = main_parser.parse_args()
    logger = logging.getLogger('preprocess.' + args.command)

    command_handlers.get(args.command)(args)