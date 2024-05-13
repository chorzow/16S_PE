### test of OTU table
import sys
import os

PWD = sys.argv[1]

with open(os.path.join(PWD, 'OTUs_output/OTU_table.tsv'), 'r') as my_file:
    with open(os.path.join(PWD, 'results/OTU_table_filtered.tsv'), 'w') as new_file:
        first_line = my_file.readline()
        new_file.write(f'{first_line}')
        for line in my_file:
            seq = line.split('\t')
            OTU = seq[1]
            if len(OTU) > 140:
                new_file.write(f'{line}')
            else:
                print(OTU)
