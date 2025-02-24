# BaseJumper BJ-Expression

The BJ-Expression pipeline is a scalable and reproducible bioinformatics pipeline to process RNAseq data and assess transcript-level and gene-level quantification. The pipeline supports both single-end and paired-end data. The pipeline takes raw sequencing data in the form of FASTQ files and performs down-sampling (randomly selecting a fixed, smaller number of reads from the full set of reads) and adapter trimming of FASTQ files. The pipeline then performs transcript-level quantification using the pseudo-alignment method Salmon as well as gene-level quantification using STAR (Spliced Transcripts Alignment to a Reference) and HTSeq.

# Pipeline Overview

The following are the steps and tools that pipeline uses to perform the analyses:

- Subsample the paired-end reads to 200,000 reads using SEQTK SAMPLE to compare metrics across samples
- Evaluate sequencing quality using FASTP and trim/clip reads
- Perform transcript-level quantification using the pseudo-alignment method implemented in SALMON
- Perform splice-aware alignment using STAR
- Extract primary aligned reads from STAR-based bam using SAMTOOLS
- Perform gene-level quantification from STAR alignment using the HTSEQ
- Evaluate STAR alignment (BAM) quality control using QUALIMAP
- Evaluate cell typing, custom metrics, and perform PCA using custom tools
- Aggregate the metrics across biosamples and tools to create overall pipeline statistics summary using MULTIQC

# Running Locally

Following are instructions for running BJ-Expression in a local Ubuntu server

## Install Java

```
sudo apt-get install default-jdk

java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Resources Required

For running the pipeline, a typical dataset requires 8 CPU cores and 50 GB of memory. For larger datasets, you may need to increase the resources to 16 CPU cores and 60 GB of memory. You can specify these resources in the command as follows:
```
--max_cpus 8 --max_memory 50.GB
```

## Test Pipeline Execution

All pipeline resources are publically available at `s3://bioskryb-public-data/pipeline_resources` users need not have to download this, and will be downloaded during nextflow run.

**Command**

example-

** csv input **

```
git clone https://github.com/BioSkryb/bj-expression.git
cd bj-expression
nextflow run main.nf --input_csv $PWD/tests/data/input/input.csv --max_cpus 8 --max_memory 50.GB
```

**Input Options**

The input for the pipeline is passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 4 columns: `biosampleName`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,read1,read2
Expression-test1,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R1_001.fastq.gz,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R2_001.fastq.gz
```

**Optional Modules**


This pipeline includes several optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_subsampling`: Set this to `true` to exclude the subsampling module. By default, it is set to `true`.
- `--skip_fastq_merge`: Set this to `true` to exclude the fastq merge module. By default, it is set to `false`.

**Outputs**

The pipeline saves its output files in the designated "publish_dir" directory. The bam files after htseq alignment are stored in the "secondary_analyses/alignment_htseq/" subdirectory and the metrics files are saved in the "secondary_analyses/secondary_metrics/" subdirectory. For details: [BJ-Expression Outputs](https://docs.basejumper.bioskryb.com/pipelines/secondary/bj-dna-qc/1.9.1/docs/#output-directories)

**command options**

```
  Usage:
      nextflow run main.nf [options]
  Script Options: see nextflow.config
  
    [required]
    --reads_csv         FILE    Path to input csv file

    --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                DEFAULT: GRCh38                            
                                
                                
    [optional]
    
    --publish_dir       DIR     Path to run output directory
                                DEFAULT: 
   
    --n_reads           VAL     Number of reads to sample for analysis
                                DEFAULT: 100000

    --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                DEFAULT: 75

    --min_reads         VAL       Minimum number of reads required for analysis.
                                DEFAULT: 5000
                                
    --skip_subsampling  STR     Skip Qualimap module
                                DEFAULT: false
                                
    --help              BOOL    Display help message
```

**Tool versions**

- `fastp: 0.20.1`
- `Seqtk: 1.3-r106`
- `Salmon: 1.6.0`
- `STAR: 2.7.6a`
- `QualiMap: 2.2.2-dev`
- `Samtools: 1.10`
- `HTSeq: 0.13.5`

**nf-test**

The BioSkryb BJ-Expression nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
sudo mv nf-test /usr/local/bin/
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("scrna-seq test") {

        when {
            params {
                // define parameters here. Example: 
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                {assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()},

                // Check for a match in the pipeline_metrics_summary csv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/pipeline_metrics_summary.csv")).match("pipeline_metrics_summary")},

                // Check for a match in the pipeline_metrics_summary_percents csv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/pipeline_metrics_summary_percents.csv")).match("pipeline_metrics_summary_percents")},

                // Check for a match in the df_dynamicrange_expression tsv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/df_dynamicrange_expression.tsv")).match("df_dynamicrange_expression")},

                // Verify existence of the df_gene_counts_salmon file
                {assert new File("${outputDir}/results_test/secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv").exists()},

                // Verify existence of the df_gene_counts_starhtseq file
                {assert new File("${outputDir}/results_test/secondary_analyses/quantification_htseq/df_gene_counts_starhtseq.tsv").exists()},

                // Verify existence of the bam file
                {assert new File("${outputDir}/results_test/secondary_analyses/alignment_htseq/Expression-test1.bam.bai").exists()}

            )
        }

    }
```
# Need Help?
If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References
For more information, you can refer to the following publications:

- Chung, C., Yang, X., Hevner, R. F., Kennedy, K., Vong, K. I., Liu, Y., Patel, A., Nedunuri, R., 
Barton, S. T., Noel, G., Barrows, C., Stanley, V., Mittal, S., Breuss, M. W., Schlachetzki,
J. C. M., Kingsmore, S. F., & Gleeson, J. G. (2024). Cell-type-resolved mosaicism 
reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392. [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V, Shao, D. D., 
Walsh, C. A., Yang, H. W., Johnson, M. D., & Park, P. J. (2024).
High-resolution 
detection of copy number alterations in single cells with HiScanner. BioRxiv, 
2024.04.26.587806. [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

- Zawistowski, J. S., Salas-González, I., Morozova, T. V, Blackinton, J. G., Tate, T., Arvapalli, 
D., Velivela, S., Harton, G. L., Marks, J. R., Hwang, E. S., Weigman, V. J., & West, J. A. 
A. (n.d.). Unifying genomics and transcriptomics in single cells with ResolveOME amplification
chemistry to illuminate oncogenic and drug resistance mechanisms. [https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full](https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of the standard 
quality control processes implemented through ResolveServices<sup>SM</sup>. While these 
pipelines may not be explicitly cited, they are integral to the methodologies 
described.
