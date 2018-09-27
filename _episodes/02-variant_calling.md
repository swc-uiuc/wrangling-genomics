---
title: "Variant Calling Workflow"
teaching: 35
exercises: 25
questions:
- "How do I find sequence variants between my samples and a reference genome?"
objectives:
- "Describe the steps involved in variant calling."
- "Describe the types of data formats encountered during variant calling."
- "Use command line tools to perform variant calling."
keypoints:
- "Bioinformatics command line tools are collections of commands that can be used to carry out bioinformatics analyses."
- "To use most powerful bioinformatics tools, you'll need to use the command line."
- "There are many different file formats for storing genomics data. It's important to understand these file formats and know how to convert among them."
---

# Alignment to a reference genome

We have already trimmed our reads so now the next step is alignment of our quality reads to the reference genome.

![workflow_align](../img/variant_calling_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be
using the [Burrows Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent
sequences against a large reference genome. The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

# Setting up

First we will link in the reference genome data into our `data/` directory, using a symbolic link.

~~~
$ cd ~/dc_workshop
$ ln -s /home/classroom/hpcbio/DC-genomics-2018/.dc_sampledata_lite/ref_genome/ data/
~~~
{: .bash}

We will also link in a set of trimmed FASTQ files to work with. These are small subsets of our real trimmed data, 
and will enable us to run our variant calling workflow quite quickly. 

~~~
$ ln -s /home/classroom/hpcbio/DC-genomics-2018/.dc_sampledata_lite/trimmed_fastq_small/ data/
~~~
{: .bash}

You will also need to create directories for the results that will be generated as part of this workflow. We can do this in a single
line of code because `mkdir` can accept multiple new directory
names as input.

~~~
$ mkdir -p results/sai results/sam results/bam results/bcf results/vcf
~~~
{: .bash}


> ## Installing Software
> 
> It's worth noting here that all of the software we are using for
> this workshop has been pre-installed on our remote computer. 
> This saves us a lot of time - installing software can be a 
> time-consuming and frustrating task - however, this does mean that
> you won't be able to walk out the door and start doing these
> analyses on your own computer. You'll need to install 
> the software first. Look at the [setup instructions](http://www.datacarpentry.org/wrangling-genomics/setup/) for more information on installing these software packages.
> 
{: .callout}

### Index the reference genome
Our first step is to index the reference genome for use by BWA. This 
helps speed up our alignment.

> ## To Index or Not To Index
> Indexing the reference only has to be run once. The only reason you would
> want to create a new index is if you are working with a different reference  
> genome or you are using a different tool for alignment.
{: .callout}

Make sure to load the appropriate module first:
~~~
$ module avail bwa  #Optional step: Displays all the modules with 'bwa' in its name
$ module load BWA/0.7.17-IGB-gcc-4.9.4
~~~
{: .bash}

~~~
$ mkdir data/bwa_index
$ bwa index -p data/bwa_index/ecoli_rel606 data/ref_genome/ecoli_rel606.fasta
~~~
{: .bash}

While the index is created, you will see output something like this:

~~~
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.52 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.59 sec
[main] Version: 0.7.15-r1140
[main] CMD: bwa index -p data/bwa_index/ecoli_rel606 data/ref_genome/ecoli_rel606.fasta
[main] Real time: 18.380 sec; CPU: 3.202 sec
~~~
{: .output}

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an 
aligner. BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence 
reads up to 100bp, while the other two are for sequences ranging from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such 
as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it 
is faster and more accurate.

Since we are working with short reads (36nt) we will be using BWA-backtrack. The general usage for BWA-backtrack is: 

**_$ bwa aln ref_genome.fasta input_file.fastq > output_file.sai_**

This will create a `.sai` file which is an intermediate file containing the suffix array indexes. 
    
Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default 
parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using 
and make sure the options you use are appropriate for your data.*

We're going to start by aligning the reads from just one of the 
samples in our dataset (`SRR098283.fastq`). Later, we'll be 
iterating this whole process on all of our sample files.

~~~
$ bwa aln data/bwa_index/ecoli_rel606 data/trimmed_fastq_small/SRR097977.fastq_trim.fastq > results/sai/SRR097977.aligned.sai
~~~
{: .bash}

You will see output that starts like this: 

~~~
[bwa_aln] 17bp reads: max_diff = 2
[bwa_aln] 38bp reads: max_diff = 3
[bwa_aln] 64bp reads: max_diff = 4
[bwa_aln] 93bp reads: max_diff = 5
[bwa_aln] 124bp reads: max_diff = 6
[bwa_aln] 157bp reads: max_diff = 7
[bwa_aln] 190bp reads: max_diff = 8
[bwa_aln] 225bp reads: max_diff = 9
[bwa_aln_core] calculate SA coordinate... 5.10 sec
[bwa_aln_core] write to the disk... 0.02 sec
~~~
{: .output}

## Alignment cleanup

![workflow_clean](../img/variant_calling_workflow_cleanup.png)

Post-alignment processing of the alignment file includes:

1. Converting output SAI alignment file to a BAM file
2. Sorting the BAM file by coordinate

### Convert the format of the alignment to SAM/BAM

The SAI file is not a standard alignment output file and will need to be converted into a SAM file before we can do any downstream
processing. 

#### SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not 
have time to go in detail of the features of the SAM format, the paper by 
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of
alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line
that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential
mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is 
displayed below with the different fields highlighted.

![sam_bam1](../img/sam_bam.png)


![sam_bam2](../img/sam_bam3.png)

First we will use the `bwa samse` command to convert the .sai file to SAM format. The usage for `bwa samse` is 

**_$ bwa samse ref_genome.fasta input_file.sai input_file.fastq > output_file.sam_**

The code in our case will look like: 

~~~
$ bwa samse data/bwa_index/ecoli_rel606 \
        results/sai/SRR097977.aligned.sai \
        data/trimmed_fastq_small/SRR097977.fastq_trim.fastq > \
        results/sam/SRR097977.aligned.sam
~~~
{: .bash}

Your output will start out something like this: 

~~~
[bwa_aln_core] convert to sequence coordinate... 0.70 sec
[bwa_aln_core] refine gapped alignments... 0.09 sec
[bwa_aln_core] print alignments... 0.37 sec
[bwa_aln_core] 262144 sequences have been processed.
~~~
{: .output}


> ## Multiple line commands
> 
> When typing a long command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
>
{: .callout}

Next we convert the SAM file to BAM format for use by downstream tools. We use the `samtools` program with the `view` command and tell this command that the input is in SAM format (`-S`) and to output BAM format (`-b`). Remember to load the samtools module first!

~~~
$ module load SAMtools/1.7-IGB-gcc-4.9.4
$ samtools view -S -b results/sam/SRR097977.aligned.sam > results/bam/SRR097977.aligned.bam
~~~
{: .bash}


### Sort BAM file by coordinates

Next we sort the BAM file using the `sort` command from `samtools`. Note that we need to specify where and what to call 
our output file with the `-o` option. 

~~~
$ samtools sort -o results/bam/SRR097977.aligned.sorted.bam results/bam/SRR097977.aligned.bam 
~~~
{: .bash}


> ## More Than One Way to . . . sort a SAM/BAM File
> By default samtools will sort by coordinates, but SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important
> to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require 
> differently sorted alignment files as input.*
{: .callout}

## Variant calling

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given position in an individual genome
or transcriptome, often referred to as a Single Nucleotide Variant (SNV). The call is usually accompanied by an estimate of 
variant frequency and some measure of confidence. Similar to other steps in this workflow, there are number of tools available for 
variant calling. In this workshop we will be using `bcftools`, but there are a few things we need to do before actually calling the variants.

![workflow](../img/variant_calling_workflow.png)

### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with samtools
[mpileup](http://samtools.sourceforge.net/mpileup.shtml):

~~~
$ samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
            results/bam/SRR097977.aligned.sorted.bam > results/bcf/SRR097977_raw.bcf
~~~
{: .bash}

~~~
[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
~~~
{: .output}

We have now generated a file with coverage information for every base. To identify variants, we now will use a different tool from the samtools suite called [bcftools](https://samtools.github.io/bcftools/bcftools.html).

### Step 2: Detect the single nucleotide polymorphisms (SNPs)

Identify SNPs using bcftools:

~~~
$ module load BCFtools/1.7-IGB-gcc-4.9.4
$ bcftools call -cv --ploidy 1 -Ob -o results/bcf/SRR097977_variants.bcf results/bcf/SRR097977_raw.bcf 
~~~
{: .bash}

Above we used several basic parameters, but it's always a good idea to look through all of the options in a 
program before running. There are two variant calling models available in BCFtools, and we are using the `-c` 
or consensus model which assumes biallelic calling rather multi-allelic calling. Since we have a haploid 
genome we shouldn't have to worry about allellic variants at all, so we're sticking with the simpler model.
We also specify that we have a haploid genome with the option `--ploidy 1`. `-v` specifies that we'd only
like to have variants reported. `-Ob` specifies that we want the output to be in bcf format, and the name of 
our output file is provided after `-o`. 

### Step 3: Filter and report the SNP variants in variant calling format (VCF)

Filter the SNPs for the final output in VCF format, using `vcfutils.pl`:

~~~
$ bcftools view results/bcf/SRR097977_variants.bcf | vcfutils.pl varFilter - > results/vcf/SRR097977_final_variants.vcf
~~~
{: .bash}

`bcftools view` converts the binary format of bcf files into human readable format (tab-delimited) for `vcfutils.pl` to perform
the filtering. Note that the output is in VCF format, which is a text format.

## Explore the VCF format:

~~~
$ less results/vcf/SRR097977_final_variants.vcf
~~~
{: .bash}

You will see the header (which describes the format), the time and date the file was
created, the version of bcftools that was used, the command line parameters used, and 
some additional information:

~~~
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=1.5+htslib-1.5
##samtoolsCommand=samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR097977.aligned.so
rted.bam
##reference=file://data/ref_genome/ecoli_rel606.fasta
##contig=<ID=NC_012967.1,length=4629812>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
.
.
.
.
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##bcftools_callVersion=1.5+htslib-1.5
##bcftools_callCommand=call --ploidy 1 -cv -Ob -o results/bcf/SRR097977_variants.bcf results/bcf/SRR097977_raw.bcf; Date=Tue Feb 20 13:30:53 2018
##bcftools_viewVersion=1.5+htslib-1.5
##bcftools_viewCommand=view results/bcf/SRR097977_variants.bcf; Date=Tue Feb 20 13:43:03 2018
~~~
{: .output}

Followed by information on each of the variations observed: 

~~~
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  results/bam/SRR097977.aligned.sorted.bam
NC_012967.1     9972    .       T       G       66.0052 .       DP=3;VDB=0.220755;SGB=-0.511536;MQSB=1;MQ0F=0;AF1=1;AC1=1;DP4=0,0,1,2;MQ=37;FQ=-999     GT:PL   1:96,0
NC_012967.1     10563   .       G       A       24.0219 .       DP=2;VDB=0.36;SGB=-0.453602;MQ0F=0;AF1=1;AC1=1;DP4=0,0,0,2;MQ=37;FQ=-999        GT:PL   1:54,0
NC_012967.1     81158   .       A       C       148.006 .       DP=8;VDB=0.250641;SGB=-0.636426;MQSB=1.01283;MQ0F=0;AF1=1;AC1=1;DP4=0,0,4,3;MQ=37;FQ=-999       GT:PL   1:178,0
NC_012967.1     216480  .       C       T       62.0051 .       DP=4;VDB=0.235765;SGB=-0.511536;MQSB=1;MQ0F=0;AF1=1;AC1=1;DP4=0,0,1,2;MQ=37;FQ=-999     GT:PL   1:92,0
NC_012967.1     247796  .       T       C       58.0051 .       DP=3;VDB=0.102722;SGB=-0.511536;MQSB=1;MQ0F=0;AF1=1;AC1=1;DP4=0,0,1,2;MQ=33;FQ=-999     GT:PL   1:88,0
NC_012967.1     281923  .       G       T       38.0055 .       DP=2;VDB=0.06;SGB=-0.453602;MQSB=1;MQ0F=0;AF1=1;AC1=1;DP4=0,0,1,1;MQ=37;FQ=-999 GT:PL   1:68,0
NC
~~~
{: .output}

This is a lot of information, so let's take some time to make sure we understand our output.

The first few columns represent the information we have about a predicted variation. 

| column | info |
| ------- | ---------- |
| CHROM | contig location where the variation occurs | 
| POS | position within the contig where the variation occurs | 
| ID | a `.` until we add annotation information | 
| REF | reference genotype (forward strand) | 
| ALT | sample genotype (forward strand) | 
| QUAL | Phred-scaled probablity that the observed variant exists at this site (higher is better) |
| FILTER | a `.` if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed | 

In an ideal world, the information in the `QUAL` column would be all we needed to filter out bad variant calls.
However, in reality we need to filter on multiple other metrics. 

The last two columns contain the genotypes and can be tricky to decode. 

| column | info |
| ------- | ---------- |
| FORMAT | lists in order the metrics presented in the final column | 
| results | lists the values associated with those metrics in order | 

For our file, the metrics presented are GT:PL:GQ. 

| metric | definition | 
| ------- | ---------- |
| GT | the genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. |
| PL | the likelihoods of the given genotypes |
| GQ | the Phred-scaled confidence for the genotype | 
| AD, DP | the depth per allele by sample and coverage |

The Broad Institute's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place
to learn more about VCF file format.

> ## Exercise
> 
> Use the `grep`, `cut`, and `less` commands you've learned to extract the `POS` and `QUAL` columns from your 
> output file (without the header lines). What is the position of the first variant to be called with a `QUAL` 
> value of less than 4?
>
>> ## Solution
>> 
>> ~~~
>> $ cut results/vcf/SRR097977_final_variants.vcf -f 6,2 | grep -v "##" | less
>> ~~~
>> {: .bash}
>> 
>> ~~~ 
>> POS     QUAL
>> 9972    66.0052
>> 10563   24.0219
>> 81158   148.006
>> 216480  62.0051
>> 247796  58.0051
>> 281923  38.0055
>> .
>> .
>> .
>> 4433347	54.005
>> 4616538	56.0051
>> ~~~
>> {: .output}
>>
>> Position 1375411 has a score of 3.01.
> {: .solution}
{: .challenge}

## Assess the alignment (visualization) - optional step

It is often instructive to look at your data in a genome browser. Visualisation will allow you to get a "feel" for 
the data, as well as detecting abnormalities and problems. Also, exploring the data in such a way may give you 
ideas for further analyses.  As such, visualization tools are useful for exploratory analysis. In this lesson we 
will describe two different tools for visualisation; a light-weight command-line based one and the Broad
Institute's Integrative Genomics Viewer (IGV) which requires
software installation and transfer of files.

In order for us to visualize the alignment files, we will need to index the BAM file using `samtools`:

~~~
$ samtools index results/bam/SRR097977.aligned.sorted.bam
~~~
{: .bash}

### Viewing with `tview`

[Samtools](http://www.htslib.org/) implements a very simple text alignment viewer based on the GNU
`ncurses` library, called `tview`. This alignment viewer works with short indels and shows [MAQ](http://maq.sourceforge.net/) consensus. 
It uses different colors to display mapping quality or base quality, subjected to users' choice. Samtools viewer is known to work with an 130 GB alignment swiftly. Due to its text interface, displaying alignments over network is also very fast.

In order to visualize our mapped reads we use `tview`, giving it the sorted bam file and the reference file: 

~~~
$ samtools tview results/bam/SRR097977.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta
~~~
{: .bash}

~~~
1         11        21        31        41        51        61        71        81        91        101       111       121
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATAC
..................................................................................................................................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,.............................
...................................,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................   ................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................................   ....................................      ,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................  ,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,     .......
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .............................  ,,,,,,,,,,,,,,,,,g,,,,,    ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ...........................T.......   ,,,,,,,,,,,,,,,,,,,,,,,c,          ......
......................... ................................   ,g,,,,,,,,,,,,,,,,,,,      ...........................
,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,       ..........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ................................T..  ..............................   ,,,,,,
...........................       ,,,,,,g,,,,,,,,,,,,,,,,,   ....................................         ,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................  ...................................        ....
....................................  ........................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
........................            .................................. .............................     ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ....................................        ..........................
...............................       ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................
...................................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ..................................
.................................... ,,,,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,        ,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ............................ ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
~~~
{: .output}

The first line of output shows the genome coordinates in our reference genome. The second line shows the reference
genome sequence. The third lines shows the consensus sequence determined from the sequence reads. A `.` indicates
a match to the reference sequence, so we can see that the consensus from our sample matches the reference in most
locations. That is good! If that wasn't the case, we should probably reconsider our choice of reference.

Below the horizontal line, we can see all of the reads in our sample aligned with the reference genome. Only 
positions where the called base differs from the reference are shown. You can use the arrow keys on your keyboard
to scroll or type `?` for a help menu. Type `Ctrl^C` to exit `tview`. 


### Viewing with IGV

[IGV](http://www.broadinstitute.org/igv/) is a stand-alone browser, which has the advantage of being installed locally and providing fast access. Web-based genome browsers, like [Ensembl](http://www.ensembl.org/index.html) or the [UCSC browser](https://genome.ucsc.edu/), are slower, but provide more functionality. They not only allow for more polished and flexible visualisation, but also provide easy access to a wealth of annotations and external data sources. This makes it straightforward to relate your data with information about repeat regions, known genes, epigenetic features or areas of cross-species conservation, to name just a few.

In order to use IGV, we will need to transfer some files to our local machine. We know how to do this with `scp`. 
Open a new tab in your terminal window and create a new folder. We'll put this folder on our Desktop for 
demonstration purposes, but in general you should avoide proliferating folders and files on your Desktop and 
instead organize files within a directory structure like we've been using in our `dc_workshop` directory.

~~~
$ mkdir ~/Desktop/files_for_igv
$ cd ~/Desktop/files_for_igv
~~~
{: .bash}

Now we will transfer our files to that new directory. Remember to replace `hpcbio##` with your temporary username. The commands to `scp` always go in the terminal window that is connected to your
local computer (not your Biocluster instance).

~~~
$ scp hpcbio##@biologin.igb.illinois.edu:~/dc_workshop/results/bam/SRR097977.aligned.sorted.bam* ~/Desktop/files_for_igv
$ scp hpcbio##@biologin.igb.illinois.edu:~/dc_workshop/data/ref_genome/ecoli_rel606.fasta ~/Desktop/files_for_igv
$ scp hpcbio##@biologin.igb.illinois.edu:~/dc_workshop/results/vcf/SRR097977_final_variants.vcf ~/Desktop/files_for_igv
~~~
{: .bash}

You will need to type in your Biocluster password each time you call `scp`. 

Next we need to open the IGV software. If you haven't done so already, you can download IGV from the [Broad Institute's software page](https://www.broadinstitute.org/software/igv/download), double-click the `.zip` file
to unzip it, and then drag the program into your Applications folder. 

1. Open IGV.
2. Load our reference genome file (`ecoli_rel606.fasta`) into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
3. Load our BAM file (`SRR097977.aligned.sorted.bam`) using the **"Load from File..."** option under the **"File"** pull-down menu. 
4.  Do the same with our VCF file (`SRR097977_final_variants.vcf`).

Your IGV browser should look like the screenshot below:

![IGV](../img/igv-screenshot.png)

There should be two tracks: one corresponding to our BAM file and the other for our VCF file. 

In the **VCF track**, each bar across the top of the plot shows the allele fraction for a single locus. The second bar shows
the genotypes for each locus in each *sample*. We only have one sample called here so we only see a single line. Dark blue = 
heterozygous, Cyan = homozygous variant, Grey = reference.  Filtered entries are transparent.

Zoom in to inspect variants you see in your filtered VCF file to become more familiar with IGV. See how quality information 
corresponds to alignment information at those loci.
Use [this website](http://software.broadinstitute.org/software/igv/AlignmentData) and the links therein to understand how IGV colors the alignments.

Now that we've run through our workflow for a single sample, we want to repeat this workflow for our other five
samples. However, we don't want to type each of these individual steps again five more times. That would be very
time consuming and error-prone, and would become impossible as we gathered more and more samples. Luckily, we
already know the tools we need to use to automate this workflow and run it on as many files as we want using a
single line of code. Those tools are: wildcards, for loops, and bash scripts. We'll use all three in the next 
lesson. 
