#!/usr/bin/env python3
"""
Generates a matrix of coverage at each variant per cell. Before using this script, first run make_cell_snv_lists.py.

Input should be the sorted bam
/home/meganle/software/Monopogen/apps/samtools view /pool001/dschaffe/variant-calling/SNU601_scATAC.possorted_bam.bam.1 | ./check_cell_coverage.py
"""
import fileinput
import numpy as np
import sys

in_snv = open("snv_list.txt", 'r')
in_cells = open("cells_list.txt", 'r')

snv_arr = []
cell_dict = {}

### convert variants to indices
# use list for SNVs because they will only be accessed sequentially
snv_line = in_snv.readline()
while(snv_line):
    snv_arr.append(snv_line.strip())
    snv_line = in_snv.readline()
in_snv.close()

# use dictionary for cells for fast index lookup
cell_idx = 0
cell_line = in_cells.readline()
while(cell_line):
    cell_dict[cell_line.strip()] = cell_idx
    cell_line = in_cells.readline()
    cell_idx += 1
in_cells.close()

cov_mat = np.zeros((len(cell_dict), len(snv_arr)))

with fileinput.input() as f_input:
    # make two passes, one for dimensions and another 
    for line in f_input:
        line_split = line.split()

        curr_var = line_split[2] + ':' + line_split[3] + ':'
        print(line, end='')

