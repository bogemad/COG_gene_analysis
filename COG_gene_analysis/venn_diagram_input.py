#!/usr/bin/env python

import sys, csv, os, itertools
from collections import defaultdict
import argparse

def write_gene_names_for_interest_groups(row):
    list = [row[0]]
    for item in row[1:]:
        if item == '':
            continue
        genes = item.split(', ')
        list += genes
    return list

def generate(orthogroups, unassigned, baseout):
    iso_d = defaultdict(int)
    out = os.path.join(baseout, "venn_output")
    os.makedirs(out, exist_ok = True)
    
    with open(orthogroups) as infile:
        samples = infile.readline().strip().split('\t')[1:]

    interest_groups = []
    for i in range(1, len(samples) + 1):
        interest_groups += itertools.combinations(samples, i)

    file_l = []
    csv_d = {}
    for i, group in enumerate(interest_groups):
        file_l.append(os.path.join(out, "{}.csv".format(".".join(group))))
        csv_d[interest_groups[i]] = csv.writer(open(file_l[i], 'w'))

    for file in (orthogroups, unassigned):
        with open(file) as infile:
            cr = csv.reader(infile, delimiter='\t')
            for i, row in enumerate(cr):
                if i == 0:
                    column_names = row[1:]
                    continue
                iso_list = []
                for j, item in enumerate(row[1:]):
                    if item.strip() == '':
                        continue
                    iso_list.append(column_names[j])
                outrow = write_gene_names_for_interest_groups(row)
                csv_d[tuple(iso_list)].writerow(outrow)
                iso_d["&".join(iso_list)] += len(write_gene_names_for_interest_groups(row)) - 1

    with open(os.path.join(out, 'theileria.combinations.txt'), 'w') as outfile:
        for key in sorted(iso_d.keys()):
            outfile.write('{}\t{}\n'.format(key, iso_d[key]))
    return file_l

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='orthogroups', help='Orthogroups.tsv file from orthofinder')
    parser.add_argument('-u', action='store', dest='unassigned', help='Orthogroups_UnassignedGenes.tsv file from orthofinder')
    parser.add_argument('-o', action='store', dest='out', help='Output directory')
    args = parser.parse_args()
    generate(args.orthogroups, args.unassigned, args.out)
