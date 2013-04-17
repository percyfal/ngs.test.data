import os
import subprocess
import logbook

LOG = logbook.Logger(__name__)

def _index_bwa(fn, label="bwa"):
    """Index bwa"""
    LOG.info("Indexing {} with bwa".format(fn))
    outdir = os.path.join(os.path.dirname(fn), os.pardir,label)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.lexists(os.path.join(outdir, os.path.basename(fn))):
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.amb".format( os.path.basename(fn)))):
        LOG.info("{} exists; not doing anything".format(fn))
        return os.path.join(outdir, os.path.basename(fn))
    cl = ["bwa", "index", os.path.abspath(os.path.join(outdir, os.path.basename(fn)))]
    subprocess.check_call(cl)
    LOG.info("Finished indexing {} with bwa".format(fn))
    return os.path.join(outdir, os.path.basename(fn))

def _index_bowtie(fn, label="bowtie"):
    """Index bowtie"""
    LOG.info("Indexing {} with bowtie".format(fn))
    outdir = os.path.join(os.path.dirname(fn), os.pardir, label)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.lexists(os.path.join(outdir, os.path.basename(fn))):
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.1.ebwt".format( os.path.splitext(os.path.basename(fn))[0]))):
        LOG.info("{} exists; not doing anything".format(fn))
        return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]
    cl = ["bowtie-build", os.path.abspath(os.path.join(outdir, os.path.basename(fn))), os.path.splitext(os.path.abspath(os.path.join(outdir, os.path.basename(fn))))[0]]
    subprocess.check_call(cl)
    LOG.info("Finished indexing {} with bowtie".format(fn))
    return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]

def _index_bowtie2(fn, label="bowtie2"):
    """Index bowtie2"""
    outdir = os.path.join(os.path.dirname(fn), os.pardir, label)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(os.path.join(outdir, os.path.basename(fn))):
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.1.bt2".format( os.path.splitext(os.path.basename(fn))[0]))):
        return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]
    cl = ["bowtie2-build", os.path.abspath(os.path.join(outdir, os.path.basename(fn))), os.path.splitext(os.path.abspath(os.path.join(outdir, os.path.basename(fn))))[0]]
    subprocess.check_call(cl)
    return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]


index_fn = {'bwa':_index_bwa, 'bowtie':_index_bowtie, 'bowtie2':_index_bowtie2}

def install_index_files(index_files):
    for k, v in index_files.iteritems():
        if not os.path.exists(os.path.dirname(v['file'])):
            os.makedirs(os.path.dirname(v['file']))
        fh = open(v['file'], "w")
        fh.write(v['data'].getvalue())
        fh.close()
