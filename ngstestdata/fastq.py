import os
import subprocess
import logbook
import yaml
from Bio import SeqIO

LOG = logbook.Logger(__name__)

with open (os.path.join(os.path.dirname(__file__), os.pardir, "config", "config.yaml")) as fh:
    config = yaml.load(fh)

def convert_bam_to_fastq(bamfile, destdir, pfx="reads"):
    """Convert bam files to fastq.  Requires picard environment variable $PICARD_HOME."""
    _bam_to_fastq(bamfile, os.path.join(destdir, pfx))
    r1 = os.path.join(destdir, "{}_1.fq".format(pfx))
    r2 = os.path.join(destdir, "{}_2.fq".format(pfx))
    files = _pair_fastq_files(r1, r2, os.path.join(destdir, "seqs"))
    return files

def _bam_to_fastq(bamfile, out_prefix):
    """Convert bam to fastq file. Outputs paired reads."""
    LOG.info("Converting bam {} to fastq with prefix {}".format(bamfile, out_prefix))
    try:
        cl = ["java", "-Xmx2g", "-XX:-UseGCOverheadLimit", "-jar", os.path.join(os.getenv("PICARD_HOME", os.curdir), "SamToFastq.jar"),
              "INPUT={}".format(bamfile), "INCLUDE_NON_PF_READS=False", "F={}_1.fq".format(out_prefix), "F2={}_2.fq".format(out_prefix), "VALIDATION_STRINGENCY=SILENT", "MAX_RECORDS_IN_RAM=5000000"]
        if not os.path.exists("{}_1.fq".format(out_prefix)):
            subprocess.check_call(cl)
    except:
        LOG.warn("Failed to run SamToFastq: {}".format(cl))
        LOG.info("This is expected since the input bamfile is truncated")
        pass

def _pair_fastq_files(r1, r2, out_prefix):
    """Pair fastq files. Unfortunately the fastq files are not
    "paired". Here we loop the files and write to outfiles only if
    there are paired reads
    """
    if os.path.exists("{}_1.fastq".format(out_prefix)):
        return
    LOG.info("reading {} and {}".format(r1, r2))
    seqs1 = _read_fastq(r1)
    seqs2 = _read_fastq(r2)
    LOG.info("Read {} sequences from read1, {} sequences from read2".format(len(seqs1), len(seqs1)))
    LOG.info("Writing fastq file 1...")
    _write_fastq("{}_1.fastq".format(out_prefix), seqs1, [x.id[0:-2] for x in seqs2])
    LOG.info("Writing fastq file 2...")
    _write_fastq("{}_2.fastq".format(out_prefix), seqs2, [x.id[0:-2] for x in seqs1])
    LOG.info("Done writing files")
    return ("{}_1.fastq".format(out_prefix), "{}_2.fastq".format(out_prefix))
    

def _read_fastq(fn, numreads=config["numreads"]):
    fh = open(fn, "rU")
    i = 0
    seqs = []
    for rec in SeqIO.parse(fh, "fastq"):
        i = i + 1
        if i % 10000==0:
            LOG.info("Read {} sequences...".format(i))
        ## For reads without a mate. SamToFastq should exclude
        ## these but apparently that doesn't happen
        if not rec.id[-2] == "/":
            LOG.warning("excluding id {}".format(rec.id))
            continue
        seqs.append(rec)
        if i >= numreads:
            return seqs

def _write_fastq(fn, seqs, ids):
    fh = open(fn, "w")
    i = 0
    for rec in seqs:
        if rec.id[-2] != "/":
            continue
        if rec.id[0:-2] in ids:
            i = i + 1
            SeqIO.write(rec, fh, "fastq")
            if i % 10000==0:
                LOG.info("Wrote {} sequences...".format(i))
    fh.close()
