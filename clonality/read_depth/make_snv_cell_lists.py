#!/usr/bin/env python3
"""
Parse coverage matrix to get list of cells/mutations of interest
"""
cell_mut = open("/pool001/meganle/scATAC_cellmut_cov5p.txt", 'r')

out_snv = open("snv_list.txt", 'w')
out_cells = open("cells_list.txt", 'w')

header = cell_mut.readline().split()
for var in header:
    out_snv.write(var + '\n')

out_snv.close()

cell_line = cell_mut.readline()
while(cell_line):
    out_cells.write(cell_line.split()[0] + '\n')
    cell_line = cell_mut.readline()

out_cells.close()
cell_mut.close()
