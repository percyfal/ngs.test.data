import os
import gzip
from Bio import SeqIO
from mako.template import Template
import logbook

LOG = logbook.Logger(__name__)


CONFIG = os.path.join(os.path.dirname(__file__), os.pardir, "config")
ARCHIVE = os.path.join(os.path.dirname(__file__), os.pardir, "data", "archive")
PRODUCTION = os.path.join(os.path.dirname(__file__), os.pardir, "data", "production")
PROJECTS = os.path.join(os.path.dirname(__file__), os.pardir, "data", "projects")
RUNINFO=Template(filename=os.path.join(CONFIG, "RunInfo.mako"))

# NOTE: following code is redundant
def install_casava_archive_files(fc_fullname, fc_id, prefix, startiter=1, nseqout=1000, **kw):
    fc_dir = os.path.join(ARCHIVE, fc_fullname)
    if not os.path.exists(fc_dir):
        os.makedirs(fc_dir)

    with open(os.path.join(fc_dir, "RunInfo.xml"), "w") as fh:
        fh.write(RUNINFO.render(**{'flowcell':fc_fullname.split("_")[3], 'fc_id':"{}_{}".format(fc_fullname.split("_")[0], fc_fullname.split("_")[3]), 'date':fc_fullname.split("_")[0], 'instrument':fc_fullname.split("_")[1]}))
    basecall_stats_dir = os.path.join(fc_dir, "Unaligned", "Basecall_Stats_{}".format(fc_id))
    if not os.path.exists(basecall_stats_dir):
        os.makedirs(basecall_stats_dir)
    for d in [os.path.join(basecall_stats_dir, x) for x in ["css", "Plots"]]:
        if not os.path.exists(d):
            os.makedirs(d)
    for row in samplesheet.split("\n"):
        if row == "":
            continue
        vals = row.split(",")
        if vals[0] == "FCID":
            header = row
            continue
        if len(vals) == 0:
            continue
        outdir = os.path.join(fc_dir, "Unaligned", "Project_{}".format(vals[5]), "Sample_{}".format(vals[2]))
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        with open(os.path.join(outdir, "SampleSheet.csv"), "w") as fh:
            LOG.info("Writing to {}".format(os.path.join(outdir, "SampleSheet.csv")))
            fh.write("{}\n".format(header))
            fh.write("{}\n".format(row))
        r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        r2 = os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        if os.path.exists(r1):
            LOG.info("{} already exists: if you want to rerun file generation remove {}".format(r1, r1))
            return 
        outf1.append(r1)
        outf2.append(r2)



def install_casava_project_files(fc_fullname, fc_id, prefix, startiter=1, nseqout=1000, **kw):
    """Install casava project files for downstream processing"""
    infile = os.path.join(CONFIG, "{}.csv".format(fc_id))
    with open(infile) as fh:
        samplesheet = fh.read()
    outf1 = []
    outf2 = []
    
    for row in samplesheet.split("\n"):
        if row == "":
            continue
        vals = row.split(",")
        if vals[0] == "FCID":
            header = row
            continue
        if len(vals) == 0:
            continue
        outdir = os.path.join(PROJECTS, vals[5].replace("__", "."), vals[2], "{}_{}".format(fc_fullname.split("_")[0], fc_fullname.split("_")[3]))
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        with open(os.path.join(outdir, "SampleSheet.csv"), "w") as fh:
            LOG.info("Writing to {}".format(os.path.join(outdir, "SampleSheet.csv")))
            fh.write("{}\n".format(header))
            fh.write("{}\n".format(row))

        r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        r2 = os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        if os.path.exists(r1):
            LOG.info("{} already exists: if you want to rerun file generation remove {}".format(r1, r1))
            return 
        outf1.append(r1)
        outf2.append(r2)

    ## Write sequences
    with open("{}_1.fastq".format(prefix), "r") as fh:
        _write_sample_fastq(fh, outf1, startiter=startiter, nseqout=nseqout)
    with open("{}_2.fastq".format(prefix), "r") as fh:
        _write_sample_fastq(fh, outf2, startiter=startiter, nseqout=nseqout)

def _write_sample_fastq(fh, outfiles, startiter=0, nseqout=1000):
    i = 0
    j = 0 - startiter
    outh = []
    for of in outfiles:
        LOG.info("Opening gzip file {}".format(of))
        oh = gzip.open(of, "w")
        outh.append(oh)
    n = len(outh)
    totseqout = n * nseqout
    LOG.info("writing {} sequences per file for {} samples, {} sequences in total".format(nseqout, n, totseqout))
    for rec in SeqIO.parse(fh, "fastq"):
        j = j + 1
        if j < 0:
            continue
        SeqIO.write(rec, outh[i%n], "fastq")
        i = i + 1
        if (i % 1000 == 0):
            LOG.info("read {} sequences from {}".format(i, fh.name))
        if i > totseqout:
            break
    [h.close() for h in outh]

    

def install_flowcell_data(fc, samplesheet):
    """Install data as defined in a samplesheet"""
    install = False
    for ss in SAMPLESHEETS[k].split("\n"):
        vals = ss.split(",")
        if vals[0]=="FCID":
            continue
        outdir = os.path.join(PRODUCTION, "{}".format(vals[5].replace("__", ".")), "{}".format(vals[2]), "{}_{}".format(FLOWCELL[k].split("_")[0],FLOWCELL[k].split("_")[-1]))
        r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        r2 = os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        LOG.info("Looking for {} and {}".format(r1, r2))
        if not os.path.exists(r1) or not os.path.exists(r2):
            install = True
            break
        if install:
            LOG.info("Installing files with run_bcbb_pipeline.py for flowcell {}".format(k))
            cl = ["run_bcbb_pipeline.py", "-s", "-g", POSTPROCESS, os.path.join(ARCHIVE, FLOWCELL[k])]
            subprocess.check_call(cl)
        else:
            LOG.info("All files present; not running run_bcbb_pipeline.py")

