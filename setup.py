"""Install ngs test data sets"""
import os
import yaml
from cStringIO import StringIO

from ngstestdata.thousandg import install_1000g_test_data
from ngstestdata.fastq import convert_bam_to_fastq
from ngstestdata.genomes import install_ucsc_genome, install_phix
from ngstestdata.variantdata import install_dbsnp_entrez, install_training_data
from ngstestdata.align import install_index_files
from ngstestdata.illumina import install_casava_project_files

from setuptools import setup, find_packages

setup(name = "ngstestdata",
      version = "0.1",
      author = "Per Unneberg",
      author_email = "punneberg@gmail.com",
      description = "download test data sets for ngs analysis",
      license = "MIT",
      packages = find_packages(),
      scripts = [],
      install_requires = [
          "PyYAML >= 3.09",
          "fabric >= 1.1.1"]
      )


## Install data after setup finishes
## Problem: currently only setting up in develop mode will work since
## otherwise data is not visible to other packages
CONFIG = os.path.join(os.path.dirname(__file__), "config")
index_files = {'sam':{'file':os.path.join(CONFIG, "tool-data", "sam_fa_indices.loc"), 'data':StringIO()},
               'bwa':{'file':os.path.join(CONFIG, "tool-data", "bwa_index.loc"), 'data':StringIO()},
               'bowtie':{'file':os.path.join(CONFIG, "tool-data", "bowtie_indices.loc"), 'data':StringIO()},
               'bowtie2':{'file':os.path.join(CONFIG, "tool-data", "bowtie2_indices.loc"), 'data':StringIO()},
               'liftOver':{'file':os.path.join(CONFIG, "tool-data", "liftOver.loc"), 'data':StringIO()}
               }


def post_setup():
    with open (os.path.join(os.path.dirname(__file__), "config", "config.yaml")) as fh:
        config = yaml.load(fh)

    # sequence data
    bamfile = install_1000g_test_data(destdir="tmp", **config['1000g'])
    fqfiles = convert_bam_to_fastq(bamfile, destdir="tmp")
    for flowcell, v in config["illumina"]["flowcells"].iteritems():
        install_casava_project_files(prefix=os.path.join("tmp", "seqs"), **v)

    # reference data
    install_phix(index_files)
    install_ucsc_genome(index_files, build="hg19", chr="chr11", start=0, end=2000000)
    install_index_files(index_files)

    # variant data
    install_dbsnp_entrez()
    (omni_out, hapmap_out, mills_out) = install_training_data()

    

post_setup()
