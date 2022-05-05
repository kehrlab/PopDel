
PopDel+: Population-wide Deletion Calling
=========================================

__Important note:__ This is an unstable prototype version of PopDel that includes the calling/genotyping of duplications and inversions.
The read pair profiles created by the stable version of PopDel are not compatible with this version of PopDel and vice versa.
For the stable version of PopDel (that currently only supports the analysis of deletions) please refer to the 'master' branch of the repository.


<b>Input</b>: BAM/CRAM files (tested for up to 50,000) from short read paired end whole genome sequencing data

<b>Output</b>: Called deletions, dandem duplications and inversions in VCF file

<b>Note:</b> The default reference genome is GRCh38 (Genome Reference Consortium Human Build 38). Other human reference builds can be specified in the options. See [Specifying the reference genome](https://github.com/kehrlab/PopDel/wiki/03.-Create-profiles-with-popdel-profile#specifying-the-reference-genome). For other diploid organism or custom reference builds, it is necessary to specify user-defined sampling intervals. See [Sampling intervals for parameter estimation](https://github.com/kehrlab/PopDel/wiki/03.-Create-profiles-with-popdel-profile#sampling-options-for-parameter-estimation).  

## Quickstart
For more detailed information see the [Wiki](https://github.com/kehrlab/PopDel/wiki).

### Installation
```
git clone https://github.com/kehrlab/PopDel.git
cd PopDel
sudo make install
```
or with conda:    

`conda install -c bioconda popdel`    

<b>Note:</b> PopDel takes significantly more time for calling variants when installed via conda.

#### Step 1: Create profile
Create insert size profiles for each individual sample
```
# Create a profile for each BAM-file
popdel profile myBam1.bam
popdel profile myBam2.bam
popdel profile myBamN.bam
```
For more options see Wiki: [PopDel Profile](https://github.com/kehrlab/PopDel/wiki/03.-Create-profiles-with-popdel-profile)

#### Step 2: Call deletions, tandem duplications and inversions
Joint calling on list of all profiles
```
# Create a list of all profiles
realpath myBam*.profile > myProfiles.txt
# Run calling on all profiles
popdel call myProfiles.txt
```
For more options see Wiki: [PopDel Call](https://github.com/kehrlab/PopDel/wiki/04.-Call-deletions-with-popdel-call)

See wiki for more information on how to view the profile with [PopDel View](https://github.com/kehrlab/PopDel/wiki/06.-Inspect-profiles-with-popdel-view) and interpret the [output in VCF-format](https://github.com/kehrlab/PopDel/wiki/05.-Output-Format:-A-(modified)-VCF).

## Citation
Sebastian Niehus, Hákon Jónsson, Janina Schönberger, Eythór Björnsson, Doruk Beyter, Hannes P. Eggertsson, Patrick Sulem, Kári Stefánsson, Bjarni V. Halldórsson, Birte Kehr.    
_PopDel identifies medium-size deletions simultaneously in tens of thousands of genomes_. Nat Commun 12, 730 (2021).
doi: https://doi.org/10.1038/s41467-020-20850-5

## Version and License
```
    Last update: 2022-05-05
    PopDel version: 2.0.0-alpha
    SeqAn version: 2.1 (with HTSlib support added by Hannes P.Eggertsson)
    Author: Sebastian Niehus (Sebastian.Niehus[at]ukr.de)
```
PopDel is distributed under the GPL-3.0. Consult the accompanying [LICENSE](https://github.com/kehrlab/PopDel/blob/master/LICENSE) file for more details.
