#!/usr/bin/env python3
"""
Generates a matrix of coverage at each variant per cell. Before using this script, first run make_cell_snv_lists.py.

Input should be the sorted bam
/home/meganle/software/Monopogen/apps/samtools view /pool001/dschaffe/variant-calling/SNU601_scATAC.possorted_bam.bam.1 | ./check_cell_coverage.py
"""
import fileinput
import numpy as np
import sys

MAX_MISMATCH = 3

def is_covered(start_pos, cigar, check_pos):
    """ Checks whether aligned read covers a given position
    Args:
        start_pos: position at the beginning of the read
        cigar: CIGAR string corresponding to alignment
        check_pos: position to check against alignment

    Returns:
        True if check_pos covered, False if not
    """
    if cigar == '*':
        return False
    # keeps track of the beginning of the previous number
    num_start = 0
    i = 0
    while(i < len(cigar)):
        char = cigar[i]
        if char == 'M':
            # alignment match, so see if check_pos in subregion
            num_bases = int(cigar[num_start:i])

            if check_pos < (start_pos + num_bases):
                return True

            start_pos += num_bases
            num_start = i + 1
        elif char == 'D':
            # deletion, so move starting position forward
            num_bases = int(cigar[num_start:i])

            start_pos += num_bases
            num_start = i + 1
            
        elif char == 'I' or char == 'S' or char == 'H':
            # insertion or clipping, so go to next subregion
            num_start = i + 1
        else:
            if not char.isdigit():
                print(cigar)
                assert 0
        i += 1
    return False

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

snv_idx = 0
snv_split = snv_arr[0].split(':')
curr_chr = snv_split[0]
curr_pos = int(snv_split[1])

prog = 0
### iterate over .bam file
with fileinput.input() as f_input:
    for line in f_input:
        if prog % 1000000 == 0:
            print(prog)
        line_split = line.split()

        # skip line if not 50 bases before SNV of interest
        pos_diff = curr_pos - int(line_split[3])
        if (curr_chr == line_split[2] and (pos_diff >= 0 and pos_diff < 50)):
            # check for mismatch filter
            # (Monopogen's default is 3, but their code checks for less than 3 ...)
            NM_field = line_split[11].split(':')
            if NM_field[0] == 'NM' and int(NM_field[2]) < MAX_MISMATCH:
                if is_covered(int(line_split[3]), line_split[5], curr_pos):
                    barcode = line_split[18].split(':')[2]
                    if barcode in cell_dict:
                        cov_mat[cell_dict[barcode]][snv_idx] += 1
        else:
            # check if we need to go to the next SNV 
            # i.e., if on a new chromosome or position is past SNV
            if curr_chr != line_split[2] or pos_diff < 0:
                snv_idx += 1
                snv_split = snv_arr[snv_idx].split(':')
                curr_chr = snv_split[0]
                curr_pos = int(snv_split[1])
        prog += 1

### save matrix to file
np.savetxt("cxm_cov_mat.txt", cov_mat, delimiter=' ')
