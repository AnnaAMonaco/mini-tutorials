# Basic ChIP-seq analysis

## Quick start to cluster computing

The following walk though assumes the following directory layout, which is also the recommended one:
```
04_ChIP
├── data
│   ├── fastQC
│   ├── bamFiles
│   └── peakCalling
├── input
│   └── fastqs
├── logs
├── src
└── readMe.txt
```
Where your parent directory is numbered and named to easily find it, and contains the minimal subdirectories `data` for your output files, `input` for your fastqs, `logs` to store log files from submitted jobs, and `src` (for *source*) for your scripts. Other useful directories can be `img` for figures, `dataframes` for tables, etc. The presence of a `readMe.txt` file containing information of the steps taken during an analysis is fundamental for recording your "experiment", and so that others can know what you did without having to ask you (both lab mates and people reading your publications in the future). This is at all effects a *protocol* or *lab book* for data analysis.

Ultimately, choose any directory and subdirectory that works for you, as long as it is clear, complete, and tidy.


###### How do I submit my job?
MarIuX uses the in-house `MXQ` scheduler to submit scripts to the cluster for computing. Exerps from `mxqsub --help` with useful information abut scheduling with `MXQ`:
```
Job environment:
  -w, --workdir=DIRECTORY   set working directory (has to exist)
                                                       (default: current workdir)
  -o, --stdout=FILE         set file to capture stdout (default: '/dev/null')
  -e, --stderr=FILE         set file to capture stderr (default: <stdout>)
  -u, --umask=MASK          set mode to use as umask   (default: current umask)
  -p, --priority=PRIORITY   set priority               (default: 127)
```
I recommend always specifying your `-stderr` to a log file so you can check for errors if the run fails.

```
Job resource information:
  Scheduling is done based on the resources a job needs and 
  on the priority given to the job.

  -j, --processors=NUMBER  set number of processors    (default: 1)
  -m, --memory=SIZE        set amount of memory        (default: 2G)
  --tmpdir=SIZE            set size of MXQ_JOB_TMPDIR  (default: 10G)
  --gpu                    request a gpu
		[SIZE] may be suffixed with a combination of T, G and M
               to specify tebibytes, gibibytes and mebibytes.
               Defaults to mebibytes if no suffix is set.

  -t, --runtime=TIME       set runtime                 (default: 15m)

        [TIME] may be suffixed with a combination of y, w, d, h and m
               to specify years, weeks, days, hours and minutes
               Defaults to minutes if no suffix is set.
```
It is **very important** to set `-memory`, `--tmpdir`, and `-t` so that your job does not get killed after 15 minutes, or almost instantly for heavy duty computing (like most aligning). Unfortunately, until you get familiar with your sequencing depth, genome size, and tool it is very much a guessing game, so err on the side of larger and longer. But don't exaggerate or IT will hate you.

```
Job grouping:
  Grouping is done by default based on the jobs resource
  and priority information, so that jobs using the same
  amount of resources and having the same priority
  are grouped and executed in parallel.
  (see 'mxqadmin --help' for details on how to close groups)

  -g, --group-id=ID  Add job to group with group_id ID
  -n, --new-group    Add job to a new group if it does not match any running or waiting group
                     (default: Add job to a group based on automatic grouping)

  -a, --command-alias=NAME       set command alias  (default: <command>)
  -N, --group-name=NAME          set group name     (default: 'default')
  -P, --group-priority=PRIORITY  set group priority (default: 127)
```
These are not super necessary, but giving your job a name will help you locate it and how it's running either with `mxqdump`, or [here](http://afk.molgen.mpg.de/mxq/mxq/).


## QC and alignment
Basic steps to check for sequencing quality with `fastqc` and proceed to read alignment.
-Input:
  -Reads .fastq
  -Reference .fasta
-Output:
  -Unfiltered alignment file .bam

### QC with `fastqc`
Submit the sequencing quality control script `run_fastqc.sh`:
```bash
cd /path/to/parent/directory #make sure you're alway here before submitting

mxqsub --processors=4 --memory=10G -t 30m --tmpdir=10G -N 'chip' --stderr ./logs/$(date "+%Y.%m.%d-%H.%M.%S").fastqc.log ./src/run_fastqc.sh
```
###### The script you will be submitting:
```bash
#!/bin/bash
set -vxe

outDir=/path/to/outputDirectory

mkdir -p $outDir/fastQC
fastqc -o $outDir/fastQC $outDir/*fastq.gz
```
As this is a very fast process, you can also run it interactively directly in the command line:
```bash
# assuming that all directories already exist
outDir=/path/to/outputDirectory

fastqc -o $outDir/fastQC $outDir/*fastq.gz
```
**NOTE**: When running interactively, it is highly recommended to use a virtual screen using `screen` or `tmux`. For example, you can open an interactive screen with `screen -S <nameofchoice>` and detach from it using `ctrl+A+D`. This way, whatever process you submitted will continue running even if you lose connection, and you can reaccess the screen with `screen -R <nameofchoice>`. To see which screens you have open, use `screen -ls`.

The output is a html file containing all the statistics. Here are examples of [good Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and [bad Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) to compare your own to.


### Making a reference index
Now that you have established the (hopefully good) quality of your data, you can map it to your reference. However, you will not use a straight fasta file, but an index built on the fasta file. Each alignment tool builds its own index, and we will be using BWA for this analysis.

First, make a directory that will contain your index. This can be in you parent directory, or in a dedicated indices directory stored alongside your reference genomes (my preferred way).
```bash
# make the directory and make sure you are in it
mkdir -p /path/to/bwaIndex
cd /path/to/bwaIndex
```
Submit the indexing script `run_bwa_index.sh`:
```bash
mxqsub --processors=4 --memory=10G -t 1h --tmpdir=10G -N 'chip' --stderr ./logs/$(date "+%Y.%m.%d-%H.%M.%S").fastqc.log ./src/run_bwa_index.sh
```
###### The script you will be submitting:
```bash
#!/bin/bash
set -vxe

fasta=/path/to/reference/fasta
#prefix of the index [same as fasta name]
prefix="fasta.name"

bwa index $fasta -p $prefix
```
This step can also easily be run interactively on command line.

**NOTE**: if the prefix is not set to *exactly* the same as what is before ".fasta", the job will fail.

### Mapping the reads
Here starts the heavy duty work. Submit the trimming and alignment script `run_bwamem_wTrim_ChIP.sh`:
```bash
cd /path/to/parent/directory # make sure you are back in your parent directory

mxqsub --processors=10 --memory=50G -t 12h --tmpdir=50G -N 'chip' --stderr ./logs/$(date "+%Y.%m.%d-%H.%M.%S").bwamem.log ./src/run_bwamem_wTrim_ChIP.sh
```
The script loops through multiple samples, the names of which need to be put in a table (`sample-table.txt` below).

###### The script you will be submitting:
```bash
#!/bin/bash
set -vxe

#set variables
alignmentIndex=/path/to/bwaIndex
ID=/path/to/sample-table.txt
inDir=/path/to/fastqs
outDir=/path/to/outputDirectory

#set temporary directory
tempDir=$TMPDIR/$RANDOM.$sampleName
mkdir -p $tempDir

for i in `cat $ID`; do 
	trim_galore --illumina --basename "${i}" --paired $inDir/"${i}"_R1_001.fastq.gz $inDir/"${i}"_R2_001.fastq.gz -o $tempDir
	bwa mem -t 30 $alignmentIndex $tempDir/"${i}"_val_1.fq.gz $tempDir/"${i}"_val_2.fq.gz | samtools view -b - > $outDir/"${i}".bwa.mem.bam
done

#clean up
rm -r $tempDir
```
The resulting BAM files containg both mapping and non-mapping reads, and we want to filter them. For this step we use `samtools`, and set out filtering flag `-F` based on what we want. [This online tool](https://broadinstitute.github.io/picard/explain-flags.html) can help you figure out the numeric code needed. An additional step to properly remove duplicates using `gatk MarkDuplicates` (Picard) was added on [<u>21.10.2024</u>].

We filter the BAM file, sort it by genomic coordinate, and index it with the script `filter_sort_index_bam.sh`:
```bash
mxqsub --processors=10 --memory=50G -t 12h --tmpdir=50G -N 'chip' --stderr ./logs/$(date "+%Y.%m.%d-%H.%M.%S").bamFilter.log ./src/filter_sort_index_bam.sh
```
###### The script you will be submitting:
```bash
#!/bin/bash
set -vxe

ID=/path/to/sample-table.txt
outDir=/path/to/outputDirectory

for i in `cat $ID`; do 
	myTempDir=$TMPDIR/$RANDOM."${i}"
	mkdir -p $myTempDir
	samtools view -b -F 3332 -q 20 $outDir/"${i}".bwa.mem.bam | samtools sort --threads 20 -T $myTempDir -n -o $myTempDir/"${i}".nameSorted.bam -
	samtools fixmate $myTempDir/"${i}".nameSorted.bam $myTempDir/"${i}".nameSorted.matefixed.bam
	samtools sort --threads 20 -T $myTempDir -o $outDir/"${i}".filtered.sorted.bam $myTempDir/"${i}".nameSorted.matefixed.bam
	gatk MarkDuplicates -I $outDir/"${i}".filtered.sorted.bam -O $outDir/"${i}".filtered.sorted.rmDup.bam -M $outDir/"${i}".dupMetrics.txt --REMOVE_DUPLICATES true 
	samtools index $outDir/"${i}".filtered.sorted.rmDup.bam
done

#clean up
rm -r $myTempDir
```

#### Alternatively Bowtie2 is commonly used for ChIP-seq alignment [UNTESTED]
Making the `bowtie2` index:
```bash
bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>
```
Mapping the reads with `bowtie2`:
```bash
#!/bin/bash
set -vxe

# adapter sequences of the used libraries, below are examples, you need to know which were used in your run
adaptor_x='CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
adaptor_y='CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT'
read1=/path/to/read1.fastq
read2=/path/to/read2.fastq
sampleName="NameOfChoice"
alignmentIndex=/path/to/bowtieIndex
outDir=/path/to/outputDirectory
PE_NAME=$outDir/$sampleName
mkdir -p $outDir

#set temporary directory
tempDir=$TMPDIR/$RANDOM.$sampleName
mkdir -p $tempDir

# trim with skewer
skewer -m pe -o $tempDir/results -t 30 -x $adaptor_x -y $adaptor_y $read1 $read2
# align with bowtie2
bowtie2 -p 2 -q --local -x $alignmentIndex -1 $tempDir/results-trimmed-pair1.fastq -2 $tempDir/results-trimmed-pair2.fastq
-S | samtools view -bS - > $PE_NAME.bowtie2.bam

#clean up
rm -r $myTempDir
```

## Peak calling
We are now ready to call peaks using `macs2`; for more in depth information about the ins and outs of `macs2`, you can check out [this HBC training page](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). As in the suggested directory structure, I would recommend having a separate output directory for the bam files and for the peak calling.
```bash
mxqsub --processors=10 --memory=50G -t 12h --tmpdir=50G -N 'chip' --stderr ./logs/$(date "+%Y.%m.%d-%H.%M.%S").macs2.log ./src/run_macs2_callpeaks.sh
```
###### The script you will be submitting:
```bash
#!/bin/bash
set -vxe

input=/path/to/filteredSorted/input.bam
ID=/path/to/sample-table.txt
outDir=/path/to/outputDirectory

for i in `cat $ID`; do 
	macs2 callpeak -t $outDir/"${i}".filtered.sorted.bam -c $input -n "${i}".peaks -B --outdir $outDir -f BAM -g mm --call-summits
done
```
In the script, the option `-g` is for genome size, and it can be a numeric value -- such as `-g 1.3e+8` -- or one of the available shortcuts for three of the main model organisms: human (`hs`), mouse (`mm`), and fruit fly (`dm`). The `-B` option generates bedGraph files for ech sample, which can be used for visualisation in the genome. 

Each sample will generate 6 output files:
- `_peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
- `_peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
- `_summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended
- `_model.R`: an R script which you can use to produce a PDF image about the model based on your data and cross-correlation plot
- `_control_lambda.bdg`: bedGraph format for input sample
- `_treat_pileup.bdg`: bedGraph format for treatment sample

How many peaks were called in each sample?
```bash
wc -l *.narrowPeak
```
With the R script you can generate a pdf plot which shows the model `macs2` generated. 
```bash
Rscript-4.3 $rootName_model.R
```
The bedGraph files can also be converted to BigWig files with the UCSC tool `bedGraphToBigWig`, but first need to be coordinate sorted. We also need to provide the size of each chromosome, which can be generated from the fasta
```bash
# get chromosome sizes from reference fasta
samtool faidx input.fasta -i chromsizes > chrom.sizes
# make bigWig
sort -k1,1 -k2,2n unsorted.bedGraph | bedGraphToBigWig chrom.sizes out.sorted.bw
```

#### Little extra: how do I run several samples in parallel?
##### Option 1: submit a script with a `for` loop
You can use a `for` loop in the script you submit with `mxqsub`, which will iterate over each sample defined in a list of samples you provide. Keep in mind that each iteration will happen one at a time, so you will need to increase the memory and/or time requirements accordingly.

The script for `run_macs2_callpeaks.sh` modified to contain a `for` loop will look like this:
```bash
#! /usr/bin/bash
set -vxe

ID=/path/to/listOfSamples.txt
inDir=/path/to/inputDirectory
outDir=/path/to/outputDirectory

for i in `cat $ID`; do 
	macs2 callpeak -t $inDir/"${i}".rmdup.bam -n "${i}"_peaks -B --outdir $outDir/"${i}"-peaks -f BAM -g mm --call-summits
done
```

##### Option 2: run or submit using `parallel`
You can use `parallel` to submit multiple samples at once, provided you have a list of the samples you want to submit. You will first set the arguments and specify the script to submit with `xargs`, you can read more about the different `xarg` options [here](https://man7.org/linux/man-pages/man1/xargs.1.html), or as usual typing `xargs --help` in the command line. The list of samples then can be passed to `parallel`, which will run the script provided - as the name suggests - in parallel. The `parallel` command is followed by the actual script command in single quotes, where the variable argument (aka, your diffeerent samples) are represented by the positional replacement `{}`.
```bash
cat list-samples.txt | xargs -P 6 -n 1 bash src/run_macs2_callpeaks.sh
cat list-samples.txt | parallel --max-procs=12 'macs2 callpeak -t {}.rmdup.bam -c {}.input.bam -n {}_peaks -B --outdir data/{}_peaks -f BAM -g mm --call-summits'
```
Using a slightly modified script that looks like this:
```bash
#! /usr/bin/bash
set -vxe
macs2 callpeak -t input/"$1".rmdup.bam -n "$1"_peaks -B --outdir data/"$1"-peaks -f BAM -g mm --call-summits
```
**NOTE**: This does not *submit* a script in parallel, but runs them in parallel interactively. To submit it, you can write a script with the `parallel` lines listed above, and submit it as previously with `mxqsub`.
