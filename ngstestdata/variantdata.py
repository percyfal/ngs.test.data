import os
from Bio import Entrez
from ngstestdata import config
import logbook

LOG = logbook.Logger(__name__)

genomes = config["genomes"]

##############################
## Variation data
##############################
dbsnp_header = """##fileformat=VCFv4.0
##fileDate=20101103
##source=dbSNP
##dbSNP_BUILD_ID=132
##reference=GRCh37
##phasing=partial
##variationPropertyDocumentationUrl=ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf
##INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=VP,Number=1,Type=String,Description="Variation Property">
##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description="First SNP Build for RS">
##INFO=<ID=WGT,Number=1,Type=Integer,Description="Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more">
##INFO=<ID=VC,Number=1,Type=String,Description="Variation Class">
##INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
##INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
##INFO=<ID=TPA,Number=0,Type=Flag,Description="Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)">
##INFO=<ID=PMC,Number=0,Type=Flag,Description="Links exist to PubMed Central article">
##INFO=<ID=S3D,Number=0,Type=Flag,Description="Has 3D structure - SNP3D table">
##INFO=<ID=SLO,Number=0,Type=Flag,Description="Has SubmitterLinkOut - From SNP->SubSNP->Batch.link_out">
##INFO=<ID=NSF,Number=0,Type=Flag,Description="Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44">
##INFO=<ID=NSM,Number=0,Type=Flag,Description="Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42">
##INFO=<ID=NSN,Number=0,Type=Flag,Description="Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41">
##INFO=<ID=REF,Number=0,Type=Flag,Description="Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8">
##INFO=<ID=SYN,Number=0,Type=Flag,Description="Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3">
##INFO=<ID=U3,Number=0,Type=Flag,Description="In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53">
##INFO=<ID=U5,Number=0,Type=Flag,Description="In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55">
##INFO=<ID=ASS,Number=0,Type=Flag,Description="In acceptor splice site FxnCode = 73">
##INFO=<ID=DSS,Number=0,Type=Flag,Description="In donor splice-site FxnCode = 75">
##INFO=<ID=INT,Number=0,Type=Flag,Description="In Intron FxnCode = 6">
##INFO=<ID=R3,Number=0,Type=Flag,Description="In 3' gene region FxnCode = 13">
##INFO=<ID=R5,Number=0,Type=Flag,Description="In 5' gene region FxnCode = 15">
##INFO=<ID=OTH,Number=0,Type=Flag,Description="Has other snp with exactly the same set of mapped positions on NCBI refernce assembly.">
##INFO=<ID=CFL,Number=0,Type=Flag,Description="Has Assembly conflict. This is for weight 1 and 2 snp that maps to different chromosomes on different assemblies.">
##INFO=<ID=ASP,Number=0,Type=Flag,Description="Is Assembly specific. This is set if the snp only maps to one assembly">
##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
##INFO=<ID=VLD,Number=0,Type=Flag,Description="Is Validated.  This bit is set if the snp has 2+ minor allele count based on frequency or genotype data.">
##INFO=<ID=G5A,Number=0,Type=Flag,Description=">5% minor allele frequency in each and all populations">
##INFO=<ID=G5,Number=0,Type=Flag,Description=">5% minor allele frequency in 1+ populations">
##INFO=<ID=HD,Number=0,Type=Flag,Description="Marker is on high density genotyping kit (50K density or greater).  The snp may have phenotype associations present in dbGaP.">
##INFO=<ID=GNO,Number=0,Type=Flag,Description="Genotypes available. The snp has individual genotype (in SubInd table).">
##INFO=<ID=KGPilot1,Number=0,Type=Flag,Description="1000 Genome discovery(pilot1) 2009">
##INFO=<ID=KGPilot123,Number=0,Type=Flag,Description="1000 Genome discovery all pilots 2010(1,2,3)">
##INFO=<ID=KGVAL,Number=0,Type=Flag,Description="1000 Genome validated by second method">
##INFO=<ID=KGPROD,Number=0,Type=Flag,Description="1000 Genome production phase">
##INFO=<ID=PH1,Number=0,Type=Flag,Description="Phase 1 genotyped: filtered, non-redundant">
##INFO=<ID=PH2,Number=0,Type=Flag,Description="Phase 2 genotyped: filtered, non-redundant">
##INFO=<ID=PH3,Number=0,Type=Flag,Description="Phase 3 genotyped: filtered, non-redundant">
##INFO=<ID=CDA,Number=0,Type=Flag,Description="Variation is interrogated in a clinical diagnostic assay">
##INFO=<ID=LSD,Number=0,Type=Flag,Description="Submitted from a locus-specific database">
##INFO=<ID=MTP,Number=0,Type=Flag,Description="Microattribution/third-party annotation(TPA:GWAS,PAGE)">
##INFO=<ID=OM,Number=0,Type=Flag,Description="Has OMIM/OMIA">
##INFO=<ID=NOC,Number=0,Type=Flag,Description="Contig allele not present in SNP allele list. The reference sequence allele at the mapped position is not present in the SNP allele list, adjusted for orientation.">
##INFO=<ID=WTD,Number=0,Type=Flag,Description="Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory">
##INFO=<ID=NOV,Number=0,Type=Flag,Description="Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common.">
##INFO=<ID=GCF,Number=0,Type=Flag,Description="Has Genotype Conflict Same (rs, ind), different genotype.  N/N is not included.">
"""

omni="""##fileformat=VCFv4.1
##FILTER=<ID=NOT_POLY_IN_1000G,Description="Alternate allele count = 0">
##FILTER=<ID=badAssayMapping,Description="The mapping information for the SNP assay is internally inconsistent in the chip metadata">
##FILTER=<ID=dup,Description="Duplicate assay at same position with worse Gentrain Score">
##FILTER=<ID=id10,Description="Within 10 bp of an known indel">
##FILTER=<ID=id20,Description="Within 20 bp of an known indel">
##FILTER=<ID=id5,Description="Within 5 bp of an known indel">
##FILTER=<ID=id50,Description="Within 50 bp of an known indel">
##FILTER=<ID=refN,Description="Reference base is N. Assay is designed for 2 alt alleles">
##FORMAT=<ID=GC,Number=.,Type=Float,Description="Gencall Score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[/gap/birdsuite/1kg/0.928975161471502.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false enable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=CR,Number=.,Type=Float,Description="SNP Callrate">
##INFO=<ID=GentrainScore,Number=.,Type=Float,Description="Gentrain Score">
##INFO=<ID=HW,Number=.,Type=Float,Description="Hardy-Weinberg Equilibrium">
##reference=human_g1k_v37.fasta
##source=infiniumFinalReportConverterV1.0
"""

hapmap = """##fileformat=VCFv4.1
##CombineVariants="analysis_type=CombineVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta rodBind=[/broad/shptmp/0.516962905488075.ASW.vcf, /broad/shptmp/0.516962905488075.CEU.vcf, /broad/shptmp/0.516962905488075.CHB.vcf, /broad/shptmp/0.516962905488075.CHD.vcf, /broad/shptmp/0.516962905488075.GIH.vcf, /broad/shptmp/0.516962905488075.JPT.vcf, /broad/shptmp/0.516962905488075.LWK.vcf, /broad/shptmp/0.516962905488075.MEX.vcf, /broad/shptmp/0.516962905488075.MKK.vcf, /broad/shptmp/0.516962905488075.TSI.vcf, /broad/shptmp/0.516962905488075.YRI.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub genotypemergeoption=UNSORTED variantmergeoption=UNION rod_priority_list=ASW,YRI,LWK,CHD,CHB,CEU,GIH,MKK,MEX,JPT,TSI printComplexMerges=false filteredAreUncalled=false minimalVCF=false setKey=set"
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[/broad/shptmp/ebanks//0.764768180511059.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##VariantsToVCF="analysis_type=VariantsToVCF input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=[X] excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta rodBind=[/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/rawdata/genotypes_chrX_ASW_phase3.3_consensus.b36_fwd.txt] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=/humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sample=null"
##reference=Homo_sapiens_assembly18.fasta
##source=VariantsToVCF
"""
mills = """##fileformat=VCFv4.1
##CombineVariants="analysis_type=CombineVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[./indel_hg19_051711_leftAligned.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub genotypemergeoption=PRIORITIZE filteredrecordsmergetype=KEEP_IF_ANY_UNFILTERED rod_priority_list=mills printComplexMerges=false filteredAreUncalled=false minimalVCF=false setKey=set assumeIdenticalSamples=false minimumN=1 masterMerge=false mergeInfoWithMaxAC=false"
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta rodBind=[/broad/shptmp/delangel/tmp/0.101786952306615.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##LeftAlignVariants="analysis_type=LeftAlignVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_b36_both.fasta rodBind=[./indel_hg18_051711.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##ValidateVariants="analysis_type=ValidateVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[./indel_hg19_051711_leftAligned_collapsed.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub validationType=ALL doNotValidateFilteredRecords=false warnOnErrors=true"
##VariantsToVCF="analysis_type=VariantsToVCF input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_b36_both.fasta rodBind=[./indel_hg18_051711_sorted.txt] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sample=null fixRef=true"
"""

vcfheader="{}\n".format("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]))


def install_dbsnp_entrez(build="hg19"):
    """Install a subset of snps using Entrez queries"""
    LOG.info("Installing dbsnp file for {}".format(genomes[build]['species']))
    if not genomes[build].get("dbsnp", None):
        return
    conf = genomes[build]['dbsnp']
    variationdir = os.path.join(os.path.abspath(os.curdir), genomes['dir'], genomes[build]['species'], build, "variation")
    if not os.path.exists(variationdir):
        os.makedirs(variationdir)
    fn = os.path.join(variationdir, conf['vcf'])
    if not os.path.exists(fn):
        try:
            # http://www.ncbi.nlm.nih.gov/books/NBK44454/#Search.how_do_i_search_dbsnp_for_the_tot
            ## This will actually only download a subset of snps
            handle = Entrez.esearch(db="snp", retmax=conf['retmax'], term="\"{}\"".format(conf['term']))
            record = Entrez.read(handle)
            records = []
            ## For some reason the first entries are more or less empty
            start = conf['start']
            delta = conf['delta']
            for i in xrange(start, len(record['IdList']), delta):
                LOG.info("retrieving dbsnp records {} - {}".format(i, i+delta))
                h = Entrez.efetch(db="snp", id=record['IdList'][i:i+delta], rettype="flt", retmax=delta)
                lbuffer = None
                lines = []
                while True:
                    if lbuffer:
                        l = lbuffer
                        lbuffer = None
                    else:
                        l = h.readline()
                    if l.startswith("rs") or len(lines) > 20:
                        if lines:
                            lbuffer = l
                            rec = _dbsnp_line(lines, conf['buildid'])
                            if rec is None:
                                break
                            records.append(rec)
                            lines = []
                        else:
                            lines.append(l.rstrip())
                    else:
                        lines.append(l.rstrip())
        except:
            LOG.warning("Entrez query failed")
            pass
        LOG.info("Writing file {}".format(fn))
        fh = open(fn, "w")
        fh.write(dbsnp_header)
        ## Header must be tab-separated, otherwise GATK complains...
        fh.write(vcfheader)
        for rec in sorted(records, key=lambda x: int(x.split()[1])):
            fh.write(rec)
            fh.write("\n")
        fh.close()

def _dbsnp_line(lines, buildid):
    if not lines:
        return None
    line = "".join(lines)
    try:
        m = re.search("^(rs[0-9]+).*alleles=([A-Z]+)/([A-Z]+).*chr=([0-9]+).*chr-pos=([0-9]+)", line)
    except:
        LOG.warning("regexp search failed")
        raise
    if m:
        return "\t".join(["chr{}".format(m.groups()[3]), m.groups()[4], m.groups()[0], m.groups()[1], m.groups()[2], ".", ".", "dbSNPBuildID={};VP=050000020005000000000100;WGT=1".format(buildid)])
    else:
        return None

def install_training_data(build="hg19"):
    variationdir = os.path.join(os.path.abspath(os.curdir), genomes['dir'], genomes[build]['species'], build, "variation")
    omni_out = os.path.join(variationdir, "1000G_omni2.5.vcf")
    hapmap_out = os.path.join(variationdir, "hapmap_3.3.vcf")
    mills_out = os.path.join(variationdir, "Mills_Devine_2hit.indels.vcf")
    fh = open(omni_out, "w")
    fh.write(omni)
    fh.write(vcfheader)
    fh.close()
    fh = open(hapmap_out, "w")
    fh.write(hapmap)
    fh.write(vcfheader)
    fh.close()
    fh = open(mills_out, "w")
    fh.write(mills)
    fh.write(vcfheader)
    fh.close()
    return (omni_out, hapmap_out, mills_out)
