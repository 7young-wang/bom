# bom

## Using R on Gadi

1. Load the required module on shell to be able to enter the R environment. Do `module load intel-compiler/2021.10.0` and `module load R/4.3.1`. To see if desired package is ready to use, use `module avail R` to find which versions on R is available. The intel-compiler is required to understand any icc command.

2. Try to download the required packages and library. I think the R library needs to be built from scratch for each individual user. My personal path to R library is created at: `/home/552/qw6850/R/x86_64-pc-linux-gnu-library/4.3`. This means, each time when I load and use R from my terminal, I am only accessing my personal library. 

3. To install package and all its depencencies automatically, use the command `install.packages("devtools", dependencies = TRUE)`

## Filtering out CREs
*Note that here I am following the original tutorial to filter out CREs at promoter regions or overlapping with exons.*

### Useful command: finding files
To get the annotation and chr_sizes files ready, I looked into Paola's folders. Use the command `find /g/data/zk16/ -type f -name '*fimo*'` to find files from interested folders. 

### Useful command: uploading files
For the annotation file, I also uploaded the one I downloaded from the gencode website. 
`scp ./Downloads/gencode.v46.annotation.gtf qw6850@gadi.nci.org.au:/g/data/zk16/qwang/bom`

If to upload a folder, I can use the code:
`scp -r ~/Desktop/ccre_conserve_summer/human_18_restricted_cre/ qw6850@gadi.nci.org.au:/g/data/zk16/qwang/bom/human_18_restricted_cre_no_filter`

After running the R, run the following code:

```
devtools::install_github("ewonglab/BOM_package")
library(BagOfMotifs)
input_bed<-’/g/data/zk16/qwang/bom/human_18_restricted_cre.bed’
chr_sizes<-’/g/data/zk16/cc3704/mouse_enh_grammar/BOM/assemblies/hg38.chrom.sizes’
#annot<-’/g/data/zk16/qwang/bom/gencode.v46.annotation.gtf’
annot<-‘/g/data/zk16/cc3704/venvs/venv_3.7.4/lib/python3.7/site-packages/deeptoolsintervals/test/GRCh38.84.gtf.gz’

filterCREs(inputBedFile = input_bed, annotFile = annot, chrSizes = chr_sizes,
                        u = 1000, d = 1000,
                        keep_proximal = FALSE,
                        remove_proximal = TRUE,
                        non_exonic = TRUE, out_bed = "human_18_restricted_cre_filt.bed",
                        ovr_dir=TRUE,
                        celloutputDir =  "/g/data/zk16/qwang/bom/")
```

## Get the FASTA files

Again, enter the R environment and run the following code. Ideally, I should wrap up all these to a .sh file and submit the job. But since this step is relatively easy (and quick) I can just run on the terminal.

```
>BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
or
>BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
>library(BagOfMotifs)
>human <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
>generateAllFasta(bedDir = "/g/data/zk16/qwang/bom/human_18_restricted_cre_no_filter/bed_files", genome = human, fastaDir="/g/data/zk16/qwang/bom/human_18_restricted_cre_no_filter/fasta_files")
```

## FIMO scan

This is the step that is taking up quite amount of time, especially given that I have a lot of CREs to scan. Thus, I am not using the BagOfMotifs package, but want to compile it directly from the terminal (FIMO suite). If I can process things in parallel, I can speed it up.

But first, I need to successfully download the MEME suite to the cluster and be able to run it. Ideally I should just email and ask Paola where her FIMO is located at, but here I am trying from scratch.

Overall, I need to download the source code, install it, and add it to the shell. 

1. Download meme-5.5.0 source code from the website and upload it to my cluster folder by:
`scp ./Downloads/meme-5.5.0.tar qw6850@gadi.nci.org.au:/home/552/qw6850`. Since it is a tar file, unzip it by `tar -xvf meme-5.5.0.tar`

2. Install it following the MEME website
```
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt 
make 
make test 
make install
```

During this step, I have encountered the `error: Python interpreter is too old` error, which is solved by loading the `python3-as-python` module to the shell environment. I have also encountered the automake problem, but solved it by following the guideline provided with the error message.

When MEME is successfully installed, I can see the folder in my HOME directory. When I open it, I can see the `meme/bin`.

3. To add `meme/bin` to the shell path: in `~/privatemodules` folder, write(nano) a modulefile for meme. The file is named meme.
```
#%Module1.0
prepend-path PATH /home/552/qw6850/meme/bin
```

Then, `module load use.own`

7. To load the software: `module load meme`

8. To start the parallel FIMO scan, write a shell script:
```
#!/bin/bash

#PBS -l storage=scratch/zk16+gdata/zk16
#PBS -l wd=/g/data/zk16/qwang/bom
#PBS -M qwang.88@berkeley.edu
#PBS -m ae
#PBS -N FIMO_less
#PBS -e /g/data/zk16/qwang/bom/log
#PBS -o /g/data/zk16/qwang/bom/log
#PBS -l ncpus=96
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -P zk16
#PBS -l walltime=10:00:00

module load use.own
module load meme
module load parallel/20191022

cd /g/data/zk16/qwang/bom/human_18_restricted_cre_no_filter

ls *.fa | parallel fimo --thresh 0.001 --o {.} /g/data/zk16/cc3704/tools/R/bom_tests/BOM_package/inst/extdata/gimme.vertebrate.v5.0.meme {}
```

9. Submit `qsub fimo_scan.sh` and check its stats by `qstat -swx 120592068`.

## Organize the files

Useful commands:
- `rm -r directory` and `rm file`
- `mkdir folder`
- `mv *.bed bed_files/` move everything ending .bed in the current directory to the directory named `bed_files`

Now, in my bom folder, I have everything organized into the original bed files (one cell type each), and the fasta files generated from these beds, and the fimo scanning results. 

## Multimodel Training

Adapted the `make_multinomial_tab.R` script from Paola. Input three arguments: the directory where all the celltype/fimo.tsv are stored, the q_value, and the output path. I will use `q_value = 0.5` as in the tutorial. As the files are quite large, I would submit it to Gadi. Though in theory parallel will be the best option, since I am still not quite fluent with it I will start by just submitting. The following script will aggregate the motif counting table, and then submit that exact table to training. Note that I have commented out the `ggplot2` library from the train_multi.R file, as it is not used and is in some kind of conflict with other libararies.

```
#!/bin/bash

#PBS -l storage=scratch/zk16+gdata/zk16
#PBS -l wd=/g/data/zk16/qwang/bom
#PBS -M qwang.88@berkeley.edu
#PBS -m ae
#PBS -N multiple_count
#PBS -e /g/data/zk16/qwang/bom/log
#PBS -o /g/data/zk16/qwang/bom/log
#PBS -l ncpus=96
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -P zk16
#PBS -l walltime=10:00:00

module load R/4.3.1
module load intel-compiler/2021.10.0

cd /g/data/zk16/qwang/bom

Rscript make_multinomial_tab.R "human_18_restricted_cre_no_filter/fimo_files/" 0.5 "human_18_restricted_cre_no_filter/multi_count.tsv"

Rscript train_multi.R "human_18_restricted_cre_no_filter/multi_count.tsv" "human_18_restricted_cre_no_filter/train.rds"

Rscript predict_multi.R "human_18_restricted_cre_no_filter" "multi_count.tsv" "train.rds"

```

The result is stored in `121035381.gadi-pbs.OU` file under the log folder. It took me 1:30 hr, and the training result is 

```
 "Training dataset:"
[1] 95484  1456
[1] "Validation dataset:"
[1] 31828  1456
[1] "Test dataset:"
[1] 31829  1456
[1] "Training dataset after filter:"
[1] 95484  1439

Stopping. Best iteration:
[6734]  train-mlogloss:0.000174 validation-mlogloss:0.000496
```

The result is suprisingly good, if I have done everything correctly. Out of 31829 test data, only 4 got the wrong label. 

## Multimodel training -- with existing motif counts

Once we have gathered the motif counts information, we can try it with different ways of celltype classification. For this, we only need to slightly edit `multi_count.tsv` such that we get right cell types matched.

Then, skip the motif counting step and start training and prediction right with the script

```
#!/bin/bash

#PBS -l storage=scratch/zk16+gdata/zk16
#PBS -l wd=/g/data/zk16/qwang/bom
#PBS -M qwang.88@berkeley.edu
#PBS -m ae
#PBS -N multiple_count
#PBS -e /g/data/zk16/qwang/bom/log
#PBS -o /g/data/zk16/qwang/bom/log
#PBS -l ncpus=96
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -P zk16
#PBS -l walltime=10:00:00

module load R/4.3.1
module load intel-compiler/2021.10.0

cd /g/data/zk16/qwang/bom

Rscript train_multi.R "human_26/26_celltype_multi_count.tsv" "human_26/train.rds"

Rscript predict_multi.R "human_26" "26_celltype_multi_count.tsv" "train.rds"

```
The `121971562.gadi-pbs` is the successfully submitted job. It took me 16min to train and validate. The result is as follows:

```
Stopping. Best iteration:
[930]   train-mlogloss:0.094140 validation-mlogloss:0.288178

[1] TRUE
[1] 0
[1] "Training dataset:"
[1] 53364  1456
[1] "Validation dataset:"
[1] 17788  1456
[1] "Test dataset:"
[1] 17789  1456
[1] "Training dataset after filter:"
[1] 53364  1425
```

The loss is somehow much larger than in the 18 cell type case. I am wondering whether not reordering the numeric celltype values could be a problem. When tested on the test set, 6.8% got the wrong label attached. Some specific celltypes (for example CBGRC) is largely confused from others.

## Processing Mouse 26 Celltype 

Downloaded the CRE files from [CATlas](http://catlas.org/renlab_downloads/wholemousebrain/sa2.subclassv3.final.peak.srt/). Then processed locally on jupyter notebook to remove duplicates (define dulicates as those CREs that are open in more than one celltype among the 26 celltypes we are targeting at). This filtering reduce the number of CREs from 653416 to 178552. Then, upload these bed files to the folder mouse_26. Get the FASTA sequence and do FIMO scanning.

After running through the pipeline, the training was done. The accuracy looks pretty good.

```
[10000] train-mlogloss:0.000216 validation-mlogloss:0.000956
[1] TRUE
[1] 0
[1] "Training dataset:"
[1] 107131   1490
[1] "Validation dataset:"
[1] 35710  1490
[1] "Test dataset:"
[1] 35711  1490
[1] "Training dataset after filter:"
[1] 107131   1483
```


