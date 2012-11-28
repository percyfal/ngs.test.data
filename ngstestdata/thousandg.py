import os
import sys
import re
import subprocess
import logging

LOG = logging.Logger(__name__)

def install_1000g_test_data(individual, bam_ext, destdir, curlfilesize, **kw):
    """Download 1000 genomes exome data in bam format

    See ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/sequence_indices/20120522.sequence.index
    for an index of recent sequencing runs. Load data into R

    df <- read.table("20120522.sequence.index", header=TRUE, sep="\t", fill=TRUE, as.is=TRUE)
    df$INDIVIDUAL =gsub("/sequence.*", "", gsub("data/", "", df$FASTQ_FILE))
    
    and select individual based on sequencing platform

    tapply(df$INSTRUMENT_MODEL, df$INDIVIDUAL, function(x) {levels(as.factor(x))})

    Here sequencing data from individual NA21137 has (arbitrarily) been chosen for download. Sequencing was
    done at BROAD institute on a Illumina HiSeq 2000.

    Method: download the fastq file and divide sequences into batches of 10000 
    sequences to emulate different projects. Then run data through pipeline to 
    generate downstream data. Downloading partial bam files with curl is possible, but
    then SamToFastq complains about unpaired mates so in the end it is better
    to download the entire bam (500M).

    The following code extracts 200000 reads that map to the 1.2 first Mb of 
    chromosome 11.
    """
    base_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/{}".format(individual)
    bam_url = os.path.join(base_url, "exome_alignment", "{}.{}".format(individual, bam_ext))
    if not os.path.exists(destdir):
        os.mkdir(destdir)
    bamfile = os.path.join(destdir, os.path.basename(bam_url))
    smallbamfile = bamfile.replace(".bam", ".small.bam")
    if not os.path.exists(smallbamfile):
        LOG.info("downloading {} from {}".format(bamfile, base_url))
        cl = ["curl", bam_url, "-o", smallbamfile, "-r", "0-{}".format(curlfilesize)]
        subprocess.check_call(cl)
        LOG.info("finished creating {}".format(smallbamfile))
    return smallbamfile
