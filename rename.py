import os
import sys

input_dir = sys.argv[1]

files = os.listdir(os.path.join(input_dir, 'all_data'))

if all([i.endswith(('R1.fastq.gz', 'R2.fastq.gz')) for i in files]):
    print('Nothing to rename')
    sys.exit(0)

r1, r2 = {}, {}
for sample in files:
    name_split = sample.split('_')
    if name_split[-1] == '1.fastq.gz' or name_split[-2] == 'R1' or name_split[-2] == 'F':
        r1[sample] = '_'.join(name_split[:-2:]) + '_R1.fastq.gz'
    elif name_split[-1] == '2.fastq.gz' or name_split[-2] == 'R2' or name_split[-2] == 'R':
        r2[sample] = '_'.join(name_split[:-2]) + '_R2.fastq.gz'
    else:
        print(f'[RENAME] WARNING: {sample} was not classified as R1 or R2')
        sys.exit(1)

def paired():
    return len(r1) == len(r2)

if paired:
    print('Preparing dataset names...')
    rename_dict = {**r1, **r2}  # merges two dictionaries for renaming
    for source, target in rename_dict.items():
        os.rename(os.path.join(input_dir, 'all_data', source),
                  os.path.join(input_dir, 'all_data', target))
    print('[RENAME] INFO: Completed')
else:
    raise IndexError('Different number of R1 and R2 files')