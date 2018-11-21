import re

from collections import defaultdict
from multiprocessing import Pool

import numpy as np

def get_delimiter_index(x, delimiter, n):
    c = 0
    for i, char in enumerate(x):
        if char == delimiter:
            c += 1
            if c == n:
                return i
    return -1

def process_vcf_chunk_np(lines, distance, split_regex, chrom_col, pos_col, alt_col):
    positions = np.empty((len(lines),), dtype=int)
    for i, line in enumerate(lines):
        positions.put(i,
                line[get_delimiter_index(line, '\t', pos_col): get_delimiter_index(line, '\t', pos_col + 1)]
                )

    offset_positions = np.roll(positions, -1)

    result = np.abs(positions - offset_positions)

    idxs = np.nonzero(result <= distance)[0]

    cooccurances = []
    for i in idxs:
        l = lines[i]
        l2 = lines[i + 1]
        cooccuring = [
                [positions[i], l[get_delimiter_index(l, '\t', alt_col) + 1: get_delimiter_index(l, '\t', alt_col + 1)]],
                [positions[i + 1], l2[get_delimiter_index(l2, '\t', alt_col) + 1: get_delimiter_index(l2, '\t', alt_col + 1)]]
                    ]
        cooccurances.append(cooccuring)

    return cooccurances

def process_vcf_chunk(lines, distance, split_regex, chrom_col, pos_col, alt_col):
    """Process the given vcf chunk to find potential coocuring variants
    
    lines: list
        - list of strings representing vcf lines
    """
    # pull out positions
    lines = [split_regex.split(line) for line in lines]
    len_lines = len(lines)
    cooccurances = defaultdict(list)

    for i in range(0, len_lines - 1):
        x = i + 1

        cooccuring = []

        p1 = int(lines[i][pos_col])
        p2 = int(lines[x][pos_col])

        if np.abs(p2 - p1) <= distance:
            cooccuring.append([p1, lines[i][alt_col]])
            cooccuring.append([p2, lines[x][alt_col]])

            while True:
                x += 1
                if len_lines > x:
                    p2 = int(lines[x][pos_col])
                    if np.abs(p2 - p1) <= distance:
                        cooccuring.append([p2, lines[x][alt_col]])
                    else:
                        break
                else:
                    break

            cooccurances[lines[i][chrom_col]].append(cooccuring)

    return cooccurances

def job_generator(f_obj, distance, split_regex,
        chrom_col=0, pos_col=1, alt_col=4, chunk_size=1000000):
    lines = []
    for i, line in enumerate(f_obj):
        lines.append(line)
        if i % chunk_size == 0:
            yield (lines, distance, split_regex, chrom_col, pos_col, alt_col)
            # pad by last element
            lines = lines[-1:]

    # get last batch
    yield (lines, distance, split_regex, chrom_col, pos_col, alt_col)

def get_cooccuring_variants(input_vcf_fp, distance, threads=1, split_regex=re.compile(r'\t')):
    """Returns potential cooccuring variants in the input vcf that are within a given distance.

    Returns:
        - dict
            {
                <chrom>: [
                    [[p1, alt], [p2, alt], [p3, alt]],
                    [[p4, alt], [p5, alt]],
                    ...
                    ]
                ...
            }
    """

    f = open(input_vcf_fp, 'r')

    # read out header
    while True:
        line = f.readline()
        if line[:6] == '#CHROM':
            break

    with Pool(threads) as p:
        results = p.starmap(process_vcf_chunk, job_generator(f, distance, split_regex,
            chrom_col=0, pos_col=1, alt_col=4))

        d = {}
        for r in results:
            for chrom, cooccurances in r.items():
                if chrom in d:
                    d[chrom] += cooccurances
                else:
                    d[chrom] = cooccurances

        return d
