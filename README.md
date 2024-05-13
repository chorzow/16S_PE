### Instruction for 16S paired-end rRNA pipeline on CPU server

#####################################################################################
### PREREQUISITES

1. Create a parent directory for your 16S rRNA data:

	`mkdir my_folder`	# where my_folder is a name of your parent directory
	`cd my_folder`
	
2. Create a folder all_data for your *.fastq.gz files (IMPORTANT STEP):

	`mkdir all_data`	
	
3. Copy all your fastq.gz files in created `all_data` folder. **IMPORTANT:** endings of your files MUST be _R1.fastq.gz and _R2.fastq.gz.

**NB:** if your filenames don't match the pattern noted in step 3, you can use `rename.py` to prepare your filenames. **Usage:** `python3 rename.py /path/to/my_folder` 

4. Edit the config file `config.yml` for your needs. It contains the following fields:
	* `PWD`: absolute pathway to your parent directory. Example: /home/lam25/my_fancy_project/. Presence of closing slash is arbitrary.
	* `OTU_threshold`: threshold for OTU frequency filtration. Must be an integer. Usually, it is 1, but if you have more than 500 samples it is better to increase the threshold up to 10.
	* `adapters_path`: path to file with adapters for trimmomatic.
	* `trimmomatic_path`: path to trimmomatic .jar file.
	* `trim_params`: trimmomatic's parameters to use in trimming step.
	* `threads`: number of threads to use while running the pipeline.


####################
### DRY RUN

1. Activate conda base environment:

	`conda activate base`

2. Activate conda snakemake environment:

	`conda activate snakemake`

3. To check if all prerequisetes were done correct run the follow command:

	`snakemake --dry-run --snakefile Snakefile_16S_PE.smk`

4. If the dry run is completed successfully, go to the section RUN PIPELINE step 1.


####################
### RUN PIPELINE

1. Go to screen mode:

	screen (press ENTER) # more information about screen: https://linuxize.com/post/how-to-use-linux-screen/

2. Activate conda base environment:

	`conda activate base`

3. Activate conda snakemake environment:

	`conda activate snakemake`
	
4. Run pipeline 16S_snakefile_for_all.smk:
	
	`snakemake --cores 20 --rerun-incomplete --snakefile Snakefile_16S_PE.smk`
	
5. Close your screen entering CTRL+A and then D
	
6. Wait for the results (4-24 h depending on the amount of files).

7. In your parent folder the tables with OTU frequencies and taxonomy is located in "absolute_path_to_your_parent_folder"/results. These tables are used for further analysis.

8. To open your screen again write the following command and press TAB, and then ENTER:

	`screen -r` 
	
9. To kill your screen, open your screen and press CTRL+A and then press K
	
####################################################################################
### Analysis of the results
Script phyloseq.R demonstrates the example of analysis of the result tables.
More tutorials for R analysis of OTU tables:

Phyloseq: https://micca.readthedocs.io/en/latest/phyloseq.html

Alpha Diversity: https://microbiome.github.io/course_2021_radboud/alpha-diversity.html

Beta Diversity: https://microbiome.github.io/course_2021_radboud/beta-diversity.html

Differential abundance analysis: https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html
