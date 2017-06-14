import luigi
from os import path
from os import system
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from RNAseq_lib import PICARD_PATH
from RNAseq_lib import GATK_PATH
from RNAseq_lib import SNP_PLOT
from python_tools import write_obj_to_file

GATK_THREAD = '1'
GATK_NCT = '10'


class prepare(luigi.Task):
    '''
    prepare directory
    '''

    OutDir = luigi.Parameter()

    def run(self):
        log_dir = path.join(self.OutDir, 'logs')
        bam_process = path.join(self.OutDir, 'bam')

        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       bam_process])

        with self.output().open('w') as prepare_log:
            prepare_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(self.OutDir))


class add_read_groups(luigi.Task):
    '''
    add group information to bam
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir = OutDir)

    def run(self):
        each_sample_dir = path.join(OutDir, 'bam', self.sample)
        if not path.exists(each_sample_dir):
            system('mkdir -p {}/tmp'.format(each_sample_dir))

        tmp = run_cmd(['java',
                       '-jar',
                       '{}'.format(PICARD_PATH),
                       'AddOrReplaceReadGroups',
                       'I={0}/{1}.bam'.format(BamDir, self.sample),
                       'O={0}/{1}.rg_added_sorted.bam'.format(
                           each_sample_dir, self.sample),
                       'SO=coordinate',
                       'RGID={0}'.format(self.sample),
                       'RGLB=onmath',
                       'RGPL=illumina',
                       'RGPU=onmath',
                       'RGSM={0}'.format(self.sample),
                       'TMP_DIR={}/tmp'.format(each_sample_dir)])

        with self.output().open('w') as add_read_groups_log:
            add_read_groups_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.add_read_groups.log'.format(OutDir, self.sample))


class mark_duplicates(luigi.Task):
    '''
    picard markduplicated reads
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return [add_read_groups(sample=each_sample, OutDir=OutDir) for each_sample in sample_list]

    def run(self):
        mark_dup_cmd = ['java',
                        '-jar',
                        PICARD_PATH,
                        'MarkDuplicates',
                        'I={0}/bam/{1}/{1}.rg_added_sorted.bam'.format(
                            OutDir, self.sample),
                        'O={0}/bam/{1}/{1}.dedupped.bam'.format(
                            OutDir, self.sample),
                        'M={0}/bam/{1}/{1}.metrics'.format(
                            OutDir, self.sample),
                        'CREATE_INDEX=false',
                        'VALIDATION_STRINGENCY=SILENT',
                        'TMP_DIR={0}/bam/{1}/tmp'.format(OutDir, self.sample)]

        index_cmd = ['samtools',
                     'index',
                     '{0}/bam/{1}/{1}.dedupped.bam'.format(OutDir, self.sample)]

        mark_duplicates_cmd = [mark_dup_cmd, index_cmd]
        mark_duplicates_log_inf = run_cmd(mark_duplicates_cmd)

        with self.output().open('w') as mark_duplicates_log:
            mark_duplicates_log.write(mark_duplicates_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.mark_duplicates.log'.format(OutDir, self.sample))


class split_ncigar_reads(luigi.Task):
    '''
    1. splits reads into exon segments
    (getting rid of Ns but maintaining grouping information) and
    hard-clip any sequences overhanging into the intronic regions
    2. reassign mapping qualities
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return [mark_duplicates(sample=each_sample, OutDir=OutDir) for each_sample in sample_list]

    def run(self):
        split_ncigar_cmd = ['java',
                            '-Djava.io.tmpdir={0}/bam/{1}/tmp'.format(
                                OutDir, self.sample),
                            '-jar',
                            GATK_PATH,
                            '-T',
                            'SplitNCigarReads',
                            '-R',
                            '{}'.format(Ref),
                            '-I',
                            '{0}/bam/{1}/{1}.dedupped.bam'.format(
                                OutDir, self.sample),
                            '-o',
                            '{0}/bam/{1}/{1}.split.bam'.format(
                                OutDir, self.sample),
                            '-rf',
                            'ReassignOneMappingQuality',
                            '-RMQF',
                            '255',
                            '-RMQT',
                            '60',
                            '-U',
                            'ALLOW_N_CIGAR_READS',
                            '--disable_bam_indexing']

        index_cmd = ['samtools',
                     'index',
                     '{0}/bam/{1}/{1}.split.bam'.format(OutDir, self.sample)]

        split_ncigar_reads_cmd = [split_ncigar_cmd, index_cmd]
        split_ncigar_reads_log_inf = run_cmd(split_ncigar_reads_cmd)

        with self.output().open('w') as split_ncigar_reads_log:
            split_ncigar_reads_log.write(split_ncigar_reads_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.split_ncigar_reads.log'.format(OutDir, self.sample))


class snp_calling(luigi.Task):
    '''
    call snp using GATK
    '''

    OutDir = luigi.Parameter()

    def requires(self):
        return [split_ncigar_reads(sample=each_sample, OutDir=OutDir) for each_sample in sample_list]

    def run(self):
        bam_list_file = '{}/bam/bam.list'.format(OutDir)
        bam_path_list = [
            '{0}/bam/{1}/{1}.split.bam'.format(OutDir, each_sample) for each_sample in sample_list]
        with open(bam_list_file, 'w') as bam_list_file_inf:
            for each_bam in bam_path_list:
                bam_list_file_inf.write('{}\n'.format(each_bam))

        tmp = run_cmd(['java',
                       '-Djava.io.tmpdir={0}/snp/tmp'.format(OutDir),
                       '-jar',
                       '{}'.format(GATK_PATH),
                       '-T',
                       'HaplotypeCaller',
                       '-nct',
                       '{}'.format(GATK_NCT),
                       '-R',
                       '{}'.format(Ref),
                       '-I',
                       '{}'.format(bam_list_file),
                       '-o',
                       '{}/snp.raw.vcf'.format(OutDir),
                       '-dontUseSoftClippedBases',
                       '-stand_call_conf',
                       '20.0'])

        with self.output().open('w') as snp_calling_log:
            snp_calling_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/snp_calling.log'.format(OutDir))


class snp_filter(luigi.Task):
    '''
    filter snp output
    '''

    OutDir = luigi.Parameter()

    def requires(self):
        return snp_calling(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['java',
                       '-jar',
                       '{}'.format(GATK_PATH),
                       '-T',
                       'VariantFiltration',
                       '-window',
                       '35',
                       '-cluster',
                       '3',
                       '-filterName',
                       'FS',
                       '-filter',
                       'FS > 30.0',
                       '-filterName',
                       'QD',
                       '-filter',
                       'QD < 2.0',
                       '-R',
                       '{}'.format(Ref),
                       '-V',
                       '{}/snp.raw.vcf'.format(OutDir),
                       '-o',
                       '{}/snp.filter.vcf'.format(OutDir)])

        with self.output().open('w') as snp_filter_log:
            snp_filter_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/snp_filter.log'.format(OutDir))


class snp_summary(luigi.Task):

    OutDir = luigi.Parameter()

    def requires(self):

        return snp_filter(OutDir=OutDir)

    def run(self):

        snp_plot_cmds = ['bcftools stats -s - {0}/snp.filter.vcf > {0}/snp.stats'.format(OutDir)]
        snp_plot_cmds.append('plot-vcfstats -p {0}/snp_plot -s {0}/snp.stats'.format(OutDir))
        snp_plot_cmds.append('Rscript {0} --snp_stats {1}/snp_plot/tstv_by_sample.0.dat --sample_inf {2} --out_dir {1}'.format(SNP_PLOT, OutDir, SampleInf))
        snp_plot_cmd_file = path.join(OutDir, 'snp.plot.sh')
        write_obj_to_file(snp_plot_cmds, snp_plot_cmd_file)

        snp_summary_cmds = ['sh', snp_plot_cmd_file]
        snp_summary_inf = run_cmd(snp_summary_cmds)
        with self.output().open('w') as snp_summary_log:
            snp_summary_log.write(snp_summary_inf)

    def output(self):
        return luigi.LocalTarget('{}/logs/snp_summary.log'.format(OutDir))


class snp_collection(luigi.Task):

    OutDir = luigi.Parameter()
    BamDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    Ref = luigi.Parameter()

    def requires(self):
        global OutDir, BamDir, sample_list, Ref, SampleInf
        OutDir = path.abspath(self.OutDir)
        BamDir = path.abspath(self.BamDir)
        Ref = path.abspath(self.Ref)
        SampleInf = self.SampleInf
        sample_list = [each.strip().split()[1]
                       for each in open(self.SampleInf)]

        return snp_summary(OutDir=OutDir)

    def run(self):

        ignore_files = ['.ignore', 'logs', 'bam', '*.sh',
                        'snp.raw.vcf*', 'snp.stats', 'snp_plot']
        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
