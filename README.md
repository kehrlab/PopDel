PopDel - Population-wide Deletion Calling
=========================================
PopDel - Fast structural deletion calling on population-scale short-read paired-end data.

## Table of Contents
1. [Quick Start](##Quick-Start)
2. [Dependencies](##Dependencies)
3. [Installation](##Installation)
4. [Running PopDel Profile](##Running-PopDel-Profile)
5. [Running PopDel Call](##Running_PopDel-Call)
6. [Output Format](##Output-Format)
7. [Running PopDel View](##Running-PopDel-View)
8. [Parallelization](##Parallelization)
9. [Further Help](##Further-Help)
10. [Sampling Intervals for Profile Parameter Estimation](##Sampling-Intervals-for-Profile-Parameter-Estimation)
11. [Version](##Version)

## 1 - Quick Start
Install
```bash
git clone https://github.com/kehrlab/PopDel.git
cd PopDel
make all
cp popdel /usr/local/bin/.
```
Run
```bash
#Create a profile for each BAM-file
popdel profile myBam1.bam
popdel profile myBam2.bam
popdel profile myBamN.bam
#Create a list of all profiles
realpath myBam*.profile > myProfiles.txt
#Run calling on all profiles
popdel call myProfiles.txt
```

## 2 - Dependencies
PopDel has the following dependencies:
- GCC 4.9 or higher
- ZLIB

ZLIB can simply be installed using apt:
```bash
sudo apt install zlib1g-dev
```
PopDel uses SeqAn (www.seqan.de) version 2.3.1 as a header library. All required SeqAn-headers are included in the PopDel download. We advise against using your own SeqAn version for PopDel, since the header files were subject to small changes that are not (yet) part of the official SeqAn release.

## 3 - Installation
Simply run 'make' in PopDel's directory. If you want to have a debug build, run 'make debug'. For simplicity, you can copy the *popdel* executable to '/usr/local/bin/.' or any other directory in your PATH.
```bash
#in popdel-directory
make all
cp popdel /usr/local/bin/.
```

## 4 - Running PopDel Profile
The first step of PopDel is the creation of insert size profiles for each individual sample.
The sample can be given as a (single-sample) BAM-file. The file may contain multiple read groups. The file may also be compressed using gzip. The BAI-index of the BAM-File must be located in the same directory and have the same filename as the BAM-file, followed by '.bai'. The default output filename is \*BAM-FILE\*.profile. The output name and path can be changed using the option '-o'.
In the simplest case, the profiling can simply be called as:
```bash
popdel profile sample.bam
```
The profiling can also be limited to one or more regions by adding the samtools-style region(s) to the end of the program call:
```bash
popdel profile sample.bam chr21 chr1:2000000-5000000 chr4:30000000
```
This creates a single profile file for the complete chr21, chr1 from 2000000-5000000 bp and chr4 from 30000000 bp till the end of the chromosome. It is however recommended to create the profiles on the whole genomes, as you can always limit the later calling to the desired regions and this will eliminate the risk that you later have to recreate a profile if you are interested in another region.
Only the profiles are needed for the calling. Each profile should have a size of about 1-2% of the original file size. They are compressed binary files and can be viewed using the *popdel view* command (see respective [chapter](##Running-PopDel-View)).
For information on the contents of the profile also see chapter [Running PopDel View](##Running-PopDel-View).

**4.1 - Sampling options for parameter estimation**
PopDel only considers a small portion of the genome for estimating the insert size distribution. Per default it uses one well-behaved region of each chromosome (see chapter [Sampling Intervals for Profile Parameter Estimation](##Sampling-Intervals-for-Profile-Parameter-Estimation)). If the first region does not contain sufficient reads, PopDel continues sampling from the next interval until the number of reads is sufficient. The default intervals PopDel uses refer to GRCh38, for other reference genomes it might be necessary to use different sampling intervals.
A file containing these user-defined intervals can be given to PopDel by using the option '-i'. The intervals in the file have to follow the samtools-style notation, with one interval per line. For reliably results, the regions should not include regions containing abnormal sequence, like telomeric or centromeric regions. It is important, that the contig names used in the file are exactly the same as in the BAM-file (or a subset thereof). If PopDel can not sample from the default or user-defined regions, please check if the chromosomes of the BAM-file are named *chrNUM:START-END* (without leading 0's). If they don't follow this exact naming pattern, the user has to define the intervals as described above. The amount of required read-pairs for parameter estimation defaults to 50,000. This value can be modified using the option '-n'. Please note that, if the profiling is restricted to certain regions of the genome (see previous section), the sampling is only performed on the chromosomes of those regions. This behavior can be overwritten by specifying user-defined sampling intervals.

```
-i, --intervals      File with genomic intervals for parameter estimation instead of default intervals.  One closed interval per line, formatted as 'CHROM:START-END', 1-based coordinates.
-n, --min-read-num   Minimum number of read pairs for parameter estimation (per read group) Default: 50000.   
```

**4.2 - Filtering options**
The filters used by PopDel profile are preconfigured for most use cases, but can be fine-tuned by the user:
```
-d, --max-deletion-size  Maximum size of deletions. Default: 10000.
-f, --flags-set    	 Only use reads with all bits of NUM set in the bam flag. Default: 33.
-F, --flags-unset  	 Only use reads with all bits of NUM unset in the bam flag. Default: 3868.
-mq, --min-mapping-qual  Only use reads with a mapping quality of at least NUM. Default: 1.
-u, --min-unclipped      Only use reads of which at least NUM bases are not clipped. Default: 50.
-s, --min-align-score    Only use reads with an alignment score relative to read length above NUM. Default: 0.8.
```
**4.3 - Running PopDel Profile on BAM-files only containing a part of the genome**
It is possible to create profiles of BAM-files that only contain a part of the genome, e.g. only chr21.
However, there is on thing to keep in mind to avoid errors. When simply called with the default parameters, PopDel will try to sample reads from the whole genome, which will result in a very low and wrong coverage estimate since the BAM-file only contains reads for chromosome 21. To avoid this there are two possibilities.  The first is to specify the region the BAM-file contains, in this case 'chr21'. The other option is to use user defined sampling intervals that lie in adequate regions of the chromosomes in the BAM-file. If you are not interested in the coverage, you can also just ignore it, as PopDel itself will not use the value later on.

**4.4 - Examples**
Run PopDel profile on *sample.bam*, creating the ouput *sample.bam.profile*.
```bash
popdel profile sample.bam
```
Run PopDel profile on *sample.bam*, creating the output *otherfolder/othername.profile*. Use the sampling intervals given in *somewhereelse/samplingIntervals.txt*
```bash
popdel profile -i somewhereelse/samplingIntervals.txt -o otherfolder/othername.profile sample.bam
```
Run PopDel profile on the compressed *sample.bam.gz*, creating the *output sample.bam.gz.profile* and filter reads with mapping quality below 30.
```bash
popdel profile -mq 30 sample.bam.gz
```
## 5 - Running PopDel Call
After creating the profiles, *popdel call* takes a list of all profiles and performs the joint calling on all samples simultaneously. The BAM-files are no longer required for this. PopDel applies a window-wise likelihood ratio test on reads of all samples after using the empirical insert size distributions and adaptive weighting to iteratively estimate the size and frequency of the potentially present deletion(s). This approach also works for overlapping deletions or deletions with a very low allele frequency. The genotyped calls are written in VCF-format (v4.2) to the file *popdel.vcf*. The option '-o' can be uses for changing the path and the name of the output file.
In the simplest case the call is:
```bash
popdel call myProfiles.txt
```
where *myProfiles.txt* contains the paths to the previously created profiles (on file per line).
On UNIX-systems this list  can conveniently be created using the function *realpath* and redirecting the output to a file.
```bash
realpath myProfileFolder/myProfile*.profile > myProfiles.txt
```
Instead of using a file containing the profiles it is also possible to define directly as command line arguments, but only if at least two profiles are used:
```bash
popdel call myProfile1.profile myProfile2.profile [...] myProfileN.profile
```

**5.1 - General calling options**
Like the profiling, the calling can be restricted to a single or multiple regions of interest. This is done by either using the option '-r' followed by the samtools-style region (multiple times for multiple regions) or by using the option '-R' followed by the path to a file containing one region of interest per line. The profiles contain their own index, so jumping to the regions and performing the calling on them is highly efficient. Other options include:
```
-b, --buffer-size           Number of buffered windows. In range [10000..inf]. Default: 200000.
-o, --out                   Output file name. Default: popdel.vcf.
-r, --region-of-interest    Genomic region 'chr:start-end' (closed interval, 1-based index). Calling is limited to this region. Multiple regions can be defined by using the parameter -r multiple times.
-R, --ROI-file              File listing one or more regions of interest, one region per line. See parameter -r.
-n, --no-regenotyping       Outputs every potential variant window without re-genotyping.
-p, --prior-probability     Prior probability of a deletion. In range [0.0..0.9999]. Default: 0.0001.
-t, --iterations            Number of iterations in EM for length estimation. Default: 15.
-u, --unsmoothed            Disable the smoothing of the insert size histogram.
```
Note that setting a higher value for the option '-d' has no effect if the profile(s) have not been created using the same or a higher value.
The changing the value for the buffer size directly influences the memory and running time. The default value should offer a good trade-off between both resources. You can expect the required memory to roughly scale 1:1 with the window buffer, meaning an increase of the buffer by factor 10 will also increase the required memory by factor 10. In our tests this lead to a running time reduction by 1/3 thanks to less overhead. Note that the number of buffered windows should always be higher than the desired maximum deletion sizes divided by 30 to avoid splitting those variants and potentially obtaining them twice in the output.

**5.2 - Filtering options for calling**
```
-d, --max-deletion-size     Maximum size of deletions. Default: 10000.
-F, --output-failed         Also output calls which did not pass the filters.
-l, --min-init-length       Minimal deletion length at initialization of iteration. Default: standard deviation.
-m, --min-length            Minimal deletion length during iteration. Default: 95th percentile of standard deviations.
-s, --min-sample-fraction   Minimum fraction of samples which is required to have enough data in the window. In range [0..1.0]. Default: 0.75.
```

**5.4 - Examples**
Perform calling on all profiles listed in *myProfiles.txt* and write the output to *myCalls.vcf*:
```bash
popdel call -o myCalls.vcf myProfiles.txt
```
Perform calling on all profiles listed in *myProfiles.txt*, only reporting deletions between length 500 and 5000 on chr21 or chr19:45000000-55000000. Only deletions with an initialization-length of at least 450 are promoted to further iterations.
```bash
popdel call -l 450 -m 500 -d 5000 -r chr21 -r chr19:45000000-55000000 myProfiles.txt
```
## 6 - Output Format
PopDel's output is a standard VCF-4.2 file inlcuding the genotypes for every sample. Every variant is defined by its genomic position and a size estimate. The precision of the size estimate mainly depends on the 'sharpness' of the insert size distribution(s) of the samples. Because the position calculations are performed in windows of 30 bp, the the starting points of the positions have an precision of +-29 bp at best. Therefore all calls are marked as 'IMPRECISE' in the INFO field.  The *LR* value in the INFO column gives a good additional quality measure for the deletion, as the value of the QUAL field will quickly cap at 100 while the *Log-Likelihoo Ratio* has no upper limit. In fact, the QUAL value is simply a PHRED-like representation of the *LR* value.

**6.1 - Special FORMAT fields**
The *YIELD* represents how many samples could be genotyped for the deletion, regardless of carrier status.
The *Likelihood-derived Allelic Depth (LAD)* represents the number of reads that shifted the likelihood ratio in favor of the REF or ALT model (or neither).
The *Distribution-derived Allelic Depth (DAD)* is similar to the *LAD*, but bases its classification on the quantiles of the distributions. Therefore it also contains counts for read pairs that support both models, or have an insert size that is too big for the deletion model.
The *First & Last (FL)* values represent the position of the first and last read pair that acted in favor of the deletion model in the  *DAD*-calculation.
The *First to Last Distance (FLD)* is the distance between those two read pairs and should roughly correspond to the size of the deletion plus the median insert size.

**6.2 - Window-wise output**
For different applications PopDel can write the output in a window-wise fashion. This behavior is enabled by setting the flag '-n'. Consider the following deletion of length 3000:
```
chr21 2000 . N <DEL> 100 PASS IMPRECISE;SVLEN=-3000;SVTYPE=DEL;AF=0.5 [...]
```
When applying the window-wise output this becomes:
```
chr21 1970 . N <DEL> 90 PASS IMPRECISE;SVLEN=-3015;SVTYPE=DEL;AF=0.4 [...]
chr21 2000 . N <DEL> 100 PASS IMPRECISE;SVLEN=-3000;SVTYPE=DEL;AF=0.5 [...]
chr21 2030 . N <DEL> 100 PASS IMPRECISE;SVLEN=-2999;SVTYPE=DEL;AF=0.5 [...]
...
chr21 5030 . N <DEL> 90 PASS IMPRECISE;SVLEN=-3012;SVTYPE=DEL;AF=0.4 [...]
```

Note that the POS field of the file now no longer represents the starting position of the variant but that of the window. As you can see, every window that is overlapped by the deletion and that passes all tests gets reported to the output file. Further, the estimates of the deletion will slightly differ from window to window.

## 7 - Running PopDel View
If you are interested in the contents of a profile you can use *popdel view* to display them in a human-readable format. The options '-e' and '-i' allow you to inspect the header and the sampled insert size histogram(s) of all read groups in the file. The output will be written to standard output and can be piped or redirected.
```
-e, --header       Write the header.
-E, --onlyHeader   Only write the header.
-i, --histograms   Write insert size histograms.
-r, --region       Limit view to this genomic region.
 ```
 The output of a profile should look like follows:
 ```
 popdel view myProfile.profile | less
#Output:
chr21   5035777 230:-10,233:-34		4:32,85:68,103:-2,110:37
chr21   5036033 3:-37,11:0,214:15	23:110,28:-28,30:67,39:45,42:-37
chr21   5036289 4:-85,13:-96,13:-73,17:-6	89:54,93:63,106:-23,126:-43,130:12,147:12,191:106
 ```
 The first two columns together describe the starting point of window in the reference. The next column refers to the first read group of the sample, the second column to the second one and so on. The value in front of each ':' denotes the positional offset of the respective read pair from the starting point of the window, while the value after each ':' denotes its deviation from the read group's median insert size.

## 8 - Parallelization
PopDel profile relies on asynchronous I/O for reducing the time needed for profile creation. PopDel call also applies asynchronous I/O, but the process is much more computationally demanding than the profiling. The calling step itself is not parallelized. If you want to run the calling on multiple CPU's oder cluster nodes, you can easily divide the calling into batches for the different chromosomes using the option '-r' and submit each job individually.
```bash
popdel call myProfiles.txt -r chr1 -o myProfiles.chr1.vcf
popdel call myProfiles.txt -r chr2 -o myProfiles.chr2.vcf
...
popdel call myProfiles.txt -r chr22 -o myProfiles.chr22.vcf
```
Of course you can also divide the single chromosomes into multiple pieces:
```bash
popdel call myProfiles.txt -r chr1:1-10000000 -o myProfiles.chr1.1.vcf
popdel call myProfiles.txt -r chr1:10000000-20000000 -o myProfiles.chr1.2.vcf
...
```
## 9 - Further Help
Basic help:
```bash
popdel profile -h
popdel call -h
```
Advanced help:
```bash
popdel profile -H
popdel call -H
```
## 10 - Sampling Intervals for Profile Parameter Estimation
The following intervals have been selected with the target in mind that they do not contain an excess amount of irregular sequence, such as centromeric or telomeric sequence or regions of low mapability. They have been selected for the reference genome GRCh38. For the analysis of alignments to a reference genome other than GRCh38 it might be necessary to use user-defined sampling intervals (see option '-i').
```
chr1:35000000-36000000
chr2:174000000-175000000
chr3:36500000-37500000
chr4:88000000-89000000
chr5:38000000-39000000
chr6:38000000-39000000
chr7:38000000-39000000
chr8:19000000-20000000
chr9:19000000-20000000
chr10:19000000-20000000
chr11:19000000-20000000
chr12:19000000-20000000
chr13:25000000-26000000
chr14:25000000-26000000
chr15:25000000-26000000
chr16:25000000-26000000
chr17:31000000-32000000
chr18:31000000-32000000
chr19:31000000-32000000
chr20:33000000-34000000
chr21:21000000-22000000
chr22:25000000-26000000
```

## 11 - VERSION
```
    Last update: 2018-09-07
    PopDel version: 1.0
    SeqAn version: 2.3.1 (modified)
    Author: Sebastian Roskosch (Sebastian.Roskosch[at]bihealth.de)
```
