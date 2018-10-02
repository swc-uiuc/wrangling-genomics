---
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using Trimmommatic."
- "Select and set multiple options for command-line bioinformatics tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---

## Cleaning Reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We vizualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. We know that all of our samples failed at
least one of the quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to
have some reads within a sample,
or some positions (near the begining or end of reads) across all
reads that are low
quality and should be discarded. We will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic Options

Trimmomatic is a program written in the Java programming language.
You don't need to learn Java to use Trimmomatic (FastQC is also
written in Java), but the fact that it's a Java program helps
explain the syntax that is used to run Trimmomatic. First find the appropriate module,
and then load it.

~~~
$ module avail trimmomatic
$ module load Trimmomatic/0.38-Java-1.8.0_152
~~~
{: .language-bash}

Note that a message from Biocluster comes up that explains how to use
Trimmomatic. The basic command to run Trimmomatic starts like this:


**_$ java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar_**


`java` tells our computer that we're running a Java program. `-jar`
is an option specifying that we're going to specify the location of
the Java program we want to run. The Java program itself will have
a `.jar` file extension.

That's just the basic command, however. Trimmomatic has a variety of
options and parameters. We will need to specify what options we want
to use for our analysis. Here are some of the options:


| option    | meaning |
| ------- | ---------- |
| `-threads` | Specify the number of processors you want Trimmomatic to use. |
|  `SE` or `PE`   | Specify whether your reads are single or paired end. |
|  `-phred33` or `-phred64` | Specify the encoding system for your quality scores. |

In addition to these options, there are various trimming steps
available:

| step   | meaning |
| ------- | ---------- |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

We said above that a basic command for Trimmomatic looks like this:

**_$ java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE_**


However, a complete command for Trimmomatic will look something like this:

**_$ java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE -threads 3 -phred64 SRR_1056.fastq SRR_1056_trimmed.fastq ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20_**


In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `SE` | that it will be taking a single end file as input |
| `-threads 3` | to use three computing threads to run (this will speed up our run) |
| `-phred64` | that the input file uses phred-64 encoding for quality scores |
| `SRR_1056.fastq` | the input file name |
|  `SRR_1056_trimmed.fastq` | the output file to create |
| `ILLUMINACLIP:SRR_adapters.fa`| to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
|`SLIDINGWINDOW:4:20` | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |

## Running Trimmomatic

Now we will run Trimmomatic on our data. To begin, we'll need to request more CPUs. By
default, Biocluster will assign you one CPU, but you can request more with the `-n` option.
To do this, we'll need to first exit the compute node we are on now, and then log onto a new compute
node while requesting more CPUs.

~~~
$ exit  #Exits current compute node
$ srun --pty -p classroom -n 3 bash
~~~
{: .language-bash}

Navigate to your data directory and reload the Trimmomatic module:

~~~
$ cd ~/dc_workshop/data/
$ module load Trimmomatic/0.38-Java-1.8.0_152
~~~
{: .language-bash}

Then make a directory for your output to be written to:

~~~
$ mkdir trimmed_fastq
~~~
{: .language-bash}

We are going to run Trimmomatic on one of our single-end samples. We
will use a sliding window of size 4 that will remove bases if the average
phred score of that window is below 20 (like in our example above). We will also
discard any reads that do not have at least 20 bases remaining after
this trimming step.

~~~
$ java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE untrimmed_fastq/SRR098283.fastq trimmed_fastq/SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
~~~
{: .language-bash}

Your output will look something like this:

~~~
TrimmomaticSE: Started with arguments:
 untrimmed_fastq/SRR098283.fastq trimmed_fastq/SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
Automatically using 3 threads
Quality encoding detected as phred33
Input Reads: 21564058 Surviving: 17030985 (78.98%) Dropped: 4533073 (21.02%)
TrimmomaticSE: Completed successfully
~~~
{: .output}

Once you're finished, attempt the following exercise:

> ## Exercise
>
> Use the output from your Trimmomatic command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep?
>
>> ## Solution
>> 1) 21.02%
>> 2) 78.98%
> {: .solution}
{: .challenge}

You may have noticed that Trimmomatic not only automatically used 3 threads (i.e. CPUs or processors),
but it also detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output file:

~~~
$ ls trimmed_fastq/SRR098283*
~~~
{: .language-bash}

~~~
trimmed_fastq/SRR098283.fastq_trim.fastq
~~~
{: .output}

The output file is also a FASTQ file. It should be smaller than our
input file because we've removed reads. We can confirm this:

~~~
$ ls -lh *_fastq/SRR098283*
~~~
{: .language-bash}

~~~
-rw-rw-r-- 1 hpcinstru02 hpcinstru02 3.0G Feb 19 19:47 SRR098283.fastq_trim.fastq
-rw-r--r-- 1 jholmes5 hpcbio_instr 3.9G Feb 19 17:23 ../untrimmed_fastq/SRR098283.fastq
~~~
{: .output}


We've just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly!

~~~
$ for infile in `ls untrimmed_fastq`
> do
> outfile="${infile}"_trim.fastq
> java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE "untrimmed_fastq/${infile}" "trimmed_fastq/${outfile}" SLIDINGWINDOW:4:20 MINLEN:20
> done
~~~
{: .language-bash}

There are two new parts in our `for` loop. First is the line:

**_$ for infile in \`ls untrimmed_fastq\`_**

And second is the line:

**_> outfile="${infile}"_trim.fastq_**


`infile` is the first variable in our loop and takes the value
of each of the untrimmed FASTQ files. It does this by taking values from the command
"ls untrimmed_fastq". This command must be surrounded by backticks. `outfile` is the
second variable in our loop and is defined by adding `_trim.fastq` to
the end of our input file name. Use `{}` to wrap the variable so that `_trim.fastq` will
not be interpreted as part of the variable name. In addition, quoting the shell variables is
a good practice AND necessary if your variables have spaces in them. For more, check [Bash Pitfalls](http://mywiki.wooledge.org/BashPitfalls).
There are no spaces before or after the `=`.

Go ahead and run the for loop. It should take a few minutes for
Trimmomatic to run for each of our six input files. Once it's done
running, take a look at your directory contents.

~~~
$ ls trimmed_fastq/
~~~
{: .language-bash}

~~~
SRR097977.fastq_trim.fastq  SRR098027.fastq_trim.fastq  SRR098281.fastq_trim.fastq
SRR098026.fastq_trim.fastq  SRR098028.fastq_trim.fastq  SRR098283.fastq_trim.fastq
~~~
{: .output}


> ## Exercise
> Earlier we looked at the first read in our `SRR098026.fastq` file and
> saw that it was very poor quality.
>
> ~~~
> $ head -n4 SRR098026.fastq
> ~~~
> {: .language-bash}
>
> ~~~
> @SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
> +SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> !!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
> ~~~
> {: .output}
>
> After filtering out bad reads, what is the first remaining read for
> this sample? What does its quality look like?
>
>> ## Solution
>> ~~~
>> $ head -n4 trimmed_fastq/SRR098026.fastq_trim.fastq
>> ~~~
>> {: .language-bash}
>>
>> ~~~
>> @SRR098026.342 HWUSI-EAS1599_1:2:1:3:655 length=35
>> GGATNGGCCTTGTATTTATGATTCTCNGAGTCTGT
>> +SRR098026.342 HWUSI-EAS1599_1:2:1:3:655 length=35
>> BB@B!B@AACBBABCCCCBBBBBB@@!B?B<ABB@
>> ~~~
>> {: .output}
>>
>> Comparing this with our quality scale:
>>
>> ~~~
>> Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
>>                   |         |         |         |         |
>> Quality score:    0........10........20........30........40
>> ~~~
>> {: .output}
>>
>> We can see that the scores are mostly in the 30+ range. This is
>> pretty good.
> {: .solution}
{: .challenge}

We've now completed the trimming and filtering steps of our quality
control process!


> ## Bonus Exercise (Advanced)
>
> Now that we've quality controlled our samples, they should perform
> better on the quality tests run by FastQC. Go ahead and re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming.
>
>> ## Solution
>>
>> In your Biocluster terminal window do:
>>
>> ~~~
>> $ module load FastQC
>> $ mkdir ~/dc_workshop/results/fastqc_trimmed_reads
>> $ fastqc -o ~/dc_workshop/results/fastqc_trimmed_reads ~/dc_workshop/data/trimmed_fastq/*.fastq
>> ~~~
>> {: .language-bash}
>>
>> In a new tab in your local computer terminal do:
>>
>> ~~~
>> $ mkdir ~/Desktop/fastqc_html/trimmed
>> $ scp hpcbio##@biologin.igb.illinois.edu:~/dc_workshop/results/fastqc_trimmed_reads/*.html ~/Desktop/fastqc_html/trimmed
>> $ open ~/Desktop/fastqc_html/trimmed/*.html
>> ~~~
>> {: .language-bash}
>>
>> Remember to replace hpcbio## with your own temporary username.
>>
>> Before trimming, one of the sequences gave a warning and another
>> failed the per base sequence quality test. After filtering, all
>> sequences pass that test.
> {: .solution}
{: .challenge}

{% include links.md %}
