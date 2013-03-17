import os
import subprocess
import logbook
import yaml
import gzip
import cStringIO
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

from ngstestdata.align import index_fn

LOG = logbook.Logger(__name__)

with open (os.path.join(os.path.dirname(__file__), os.pardir, "config", "config.yaml")) as fh:
    config = yaml.load(fh)
genomes = config["genomes"]

def _index_seq_files(index_files, fn, build):
    outfile = index_fn['bwa'](fn, label="bwa")
    index_files['sam']['data'].write("index\t{}\t{}\n".format(build, fn))
    index_files['bwa']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    outfile = index_fn['bowtie'](fn, label="bowtie")
    index_files['bowtie']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    return index_files

def install_phix(index_files):
    LOG.info("Installing phix")
    build = "phix"
    genomedir = os.path.join(os.path.abspath(os.curdir), genomes['dir'], genomes[build]['species'], build, "seq")
    fn = os.path.join(genomedir, "phix.fa")
    if not os.path.exists(genomedir):
        LOG.info("Creating {}".format(genomedir))
        os.makedirs(genomedir)
    if not os.path.exists(fn):
        try:
            LOG.info("Opening file {}".format(fn))
            fh = open(fn, "w")
            handle = Entrez.efetch(db="nucleotide", id="9626372", rettype="fasta", retmode="text")
            rec = "".join(handle.readlines())
            fh.write(rec)
            fh.close()
        except:
            pass
    index_files = _index_seq_files(index_files, fn, build)
    return index_files

def install_seqcap_targets(build, chr):
    LOG.info("Installing target and bait files for genome build {}, chr {}".format(build, chr))
    seqcap = genomes[build]["seqcap"][chr]
    targets = cStringIO.StringIO()
    targets.write("\n".join(seqcap["header"]))
    targets.write("\n")
    targets.write("\n".join(["{0}\t{1}\t{2}\t+\tTARGET_{3}".format(chr, x.split(",")[0].lstrip(), x.split(",")[1].lstrip(), x.split(",")[0].lstrip()) for x in seqcap["targets"]]))
    targets.write("\n")
    baits = cStringIO.StringIO()
    baits.write("\n".join(seqcap["header"]))
    baits.write("\n")
    baits.write("\n".join(["{0}\t{1}\t{2}\t+\tBAIT_{3}".format(chr, x.split(",")[0].lstrip(), x.split(",")[1].lstrip(), x.split(",")[0].lstrip()) for x in seqcap["baits"]]))
    baits.write("\n")
    seqcapdir = os.path.join(os.path.abspath(os.curdir), genomes['dir'], genomes[build]['species'], build, "seqcap")
    if not os.path.exists(seqcapdir):
        LOG.info("Creating {}".format(seqcapdir))
        os.makedirs(seqcapdir)
    with open(os.path.join(seqcapdir, "{}_targets.interval_list".format(chr)), "w") as fh:
        fh.write(targets.getvalue())
    with open(os.path.join(seqcapdir, "{}_baits.interval_list".format(chr)), "w") as fh:
        fh.write(baits.getvalue())

def install_ucsc_genome(index_files, build, chr, start=None, end=None):
    LOG.info("Installing genome build {}, chr {}".format(build, chr))
    if start and end:
        LOG.info("Installing region {}-{}".format(start, end))
    ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/chromosomes/".format(build)
    url = os.path.join(ucsc, "{}.fa.gz".format(chr))
    genomedir = os.path.join(os.path.abspath(os.curdir), genomes['dir'], genomes[build]['species'], build, "seq")
    if not os.path.exists(genomedir):
        LOG.info("Creating {}".format(genomedir))
        os.makedirs(genomedir)
    try:
        LOG.info("Downloading {} from {} with curl".format(os.path.join(genomedir, os.path.basename(url)), url))
        cl = ["curl", url, "-o", os.path.join(genomedir, os.path.basename(url))]
        if not os.path.exists(os.path.join(genomedir, os.path.basename(url))):
            subprocess.check_call(cl)
    except:
        pass
    fn = os.path.join(genomedir, os.path.basename(url).replace(".gz", ""))
    if not os.path.exists(fn):
        rec = SeqIO.read(gzip.open(os.path.join(genomedir, os.path.basename(url)), "r"), "fasta")
        outh = open(fn, "w")
        SeqIO.write(SeqRecord(rec.seq[start:end], rec.id, '', ''), outh, "fasta")
        outh.close()
    index_files = _index_seq_files(index_files, fn, build)
    return index_files
