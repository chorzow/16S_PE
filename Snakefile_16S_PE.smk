import os
import sys
from snakemake.logging import logger

configfile: "config.yml"

PWD = config['PWD']
thresh = int(config['OTU_threshold'])
logger.info(f'Workdir is {PWD}')
logger.info(f'OTU threshold is {thresh}')

filenames = [i.rstrip('fastq.gz') for i in os.listdir(os.path.join(PWD, 'all_data'))]
logger.info(filenames)
#Prepare dataset basename
Sample_list = [i[:-3] for i in filenames]
logger.info(Sample_list)


#snakemake rules
ruleorder: fastqc_1 > decompress > compress > trimmomatic_r2 > trimmomatic > fastqc_2 > dada_1


rule all:
    input:
        # expand('{PWD}/results/phyloseq.png', PWD = PWD),
        expand('{PWD}/results/all_phylogeny.tsv', PWD = PWD),
        expand('{PWD}/results/OTU_table_filtered.tsv', PWD = PWD),
        expand('{PWD}/OTUs_output/OTU_table.tsv', PWD = PWD),
        expand('{PWD}/mmseq/new_DB_clu_rep.fasta', PWD = PWD),
        expand('{PWD}/seqtabnochim.fasta', PWD = PWD),
        expand('{PWD}/Fastqc_analysis/{sample}_R1_fastqc.html', PWD = PWD, sample = Sample_list),
        expand(['{PWD}/Fastqc_analysis/trimmed/{sample}.1P_fastqc.html',
                '{PWD}/Fastqc_analysis/trimmed/{sample}.2P_fastqc.html'],
                PWD = PWD, sample = Sample_list),
        expand(['{PWD}/trimmed/{sample}.1P.fastq.gz', '{PWD}/trimmed/{sample}.1U.fastq.gz',
                '{PWD}/trimmed/{sample}.2P.fastq.gz', '{PWD}/trimmed/{sample}.2U.fastq.gz'], 
                PWD = PWD, sample = Sample_list),
        expand('{PWD}/cut/{sample}_R1.fastq.gz', PWD = PWD, sample = Sample_list),
        expand('{PWD}/cut/{sample}_R2.fastq.gz', PWD = PWD, sample = Sample_list)


rule fastqc_1:
    input:
        r1 = '{PWD}/all_data/{sample}_R1.fastq.gz',
        r2 = '{PWD}/all_data/{sample}_R2.fastq.gz'
    output:
        '{PWD}/Fastqc_analysis/{sample}_R1_fastqc.html',
        '{PWD}/Fastqc_analysis/{sample}_R2_fastqc.html',
    threads: 5
    priority: 10
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        "fastqc -t 5 -o {PWD}/Fastqc_analysis/ {input.r1}; "
        "fastqc -t 5 -o {PWD}/Fastqc_analysis/ {input.r2} "
        "&& conda activate base"

rule decompress:
    input:
        check = rules.fastqc_1.output,
        gz1 = '{PWD}/all_data/{sample}_R1.fastq.gz',
        gz2 = '{PWD}/all_data/{sample}_R2.fastq.gz'
    output:
        fq1 = '{PWD}/all_data/{sample}_R1.fastq',
        fq2 = '{PWD}/all_data/{sample}_R2.fastq'
    threads: 20
    priority: 9
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        'gzip -d {input.gz1}; '
        'gzip -d {input.gz2} '
        "&& conda activate base"


rule cut_lines:
    input:
        r1 = rules.decompress.output.fq1,
        r2 = rules.decompress.output.fq2
    output:
        r1 = '{PWD}/cut/{sample}_R1.fastq',
        r2 = '{PWD}/cut/{sample}_R2.fastq',
    threads: 20
    priority: 8
    shell:
        'head -n 400000 {input.r1} > {output.r1}; '
        'head -n 400000 {input.r2} > {output.r2} '

rule compress:
    input:
        c1 = rules.cut_lines.output.r1,  # cut1
        c2 = rules.cut_lines.output.r2,  # cut2
        f1 = rules.decompress.output.fq1,  # full1
        f2 = rules.decompress.output.fq2  # full2
    output:
        c1 = '{PWD}/cut/{sample}_R1.fastq.gz',
        c2 = '{PWD}/cut/{sample}_R2.fastq.gz'
    threads: 20
    shell:
        'gzip {input.c1} ;'
        'gzip {input.c2} ;'
        'gzip {input.f1} ;'
        'gzip {input.f2}'

rule trimmomatic_r2:
    input:
        '{PWD}/cut/{sample}_R2.fastq.gz'
    output:
        '{PWD}/trimmed_r2/{sample}_R2.fastq.gz'
    threads: 20
    shell:
        'java -jar {config[trimmomatic_path]} SE -threads {threads} -phred33 {input} {output} HEADCROP:4'

rule trimmomatic:
    input:
        r1 = '{PWD}/cut/{sample}_R1.fastq.gz',
        r2 = rules.trimmomatic_r2.output
    output:
        p1 = '{PWD}/trimmed/{sample}.1P.fastq.gz',
        u1 = '{PWD}/trimmed/{sample}.1U.fastq.gz',
        p2 = '{PWD}/trimmed/{sample}.2P.fastq.gz',
        u2 = '{PWD}/trimmed/{sample}.2U.fastq.gz'
    threads: 20
    shell:
        'java -jar {config[trimmomatic_path]} PE -threads {threads} -phred33 {input.r1} {input.r2} ' \
        '{output.p1} {output.u1} {output.p2} {output.u2} {config[trim_params]}'


rule fastqc_2:
    input:
        r1 = rules.trimmomatic.output.p1,
        r2 = rules.trimmomatic.output.p2
    output:
        '{PWD}/Fastqc_analysis/trimmed/{sample}.1P_fastqc.html',
        '{PWD}/Fastqc_analysis/trimmed/{sample}.2P_fastqc.html'
    threads: 20
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        'fastqc -t {threads} -o {PWD}/Fastqc_analysis/trimmed/ {input.r1}; '
        'fastqc -t {threads} -o {PWD}/Fastqc_analysis/trimmed/ {input.r2} '
        "&& conda activate base"


rule dada_1:
    input:
        r1 = expand('{PWD}/Fastqc_analysis/trimmed/{sample}.1P_fastqc.html', PWD = PWD, sample = Sample_list),
        r2 = expand('{PWD}/Fastqc_analysis/trimmed/{sample}.2P_fastqc.html', PWD = PWD, sample = Sample_list),
    output:
        '{PWD}/seqtabnochim.fasta'
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        "Rscript DADA2_part_1.R {PWD} {thresh} &> {PWD}/dada_1.log "
        "&& conda activate base"

rule mmseq_part_1:
    input:
        rules.dada_1.output
    output:
        '{PWD}/mmseq/new_DB_clu_rep.fasta'
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        "bash mmseq_clustering.sh {PWD} &> {PWD}/mmseq_part_1.log "
        "&& conda activate base"

rule mmseq_part_2:
    input:
        rules.mmseq_part_1.output
    output:
        '{PWD}/OTUs_output/OTU_table.tsv'
    threads: 20
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        'python3 mmseq2_clustering.py {PWD} &> {PWD}/mmseq_part_2.log '
        "&& conda activate base"

rule filter_OTUs:
    input:
        rules.mmseq_part_2.output
    output:
        '{PWD}/results/OTU_table_filtered.tsv'
    threads: 20
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        'python3 filter_OTUs.py {PWD} '
        "&& conda activate base"

rule dada_2:
    input:
        rules.filter_OTUs.output
    output:
        '{PWD}/results/all_phylogeny.tsv'
    shell:
        "set +eu "
        "&& PS1=dummy "
        "&& . $(conda info --base)/etc/profile.d/conda.sh "
        "&& conda activate seurat "
        "&& echo $CONDA_PREFIX; "
        'Rscript DADA2_part_2.R {PWD} &> {PWD}/dada_2.log '
        "&& conda activate base"

# No phyloseq package on CPU-server!!!
# rule phyloseq:
#     input:
#         rules.dada_2.output
#     output:
#         '{PWD}/results/phyloseq.png'
#     shell:
#         'Rscript phyloseq.R {PWD} &> phyloseq.log'
