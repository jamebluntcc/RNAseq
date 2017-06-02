import luigi
from os import path
from os import system
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from RNAseq_lib import STAR_MAPPING_STATS
from RNAseq_lib import STAR_MAPPING_STATS_PLOT

STAR_THREAD = 8


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        mapping_dir = path.join(OutDir, 'mapping_dir')
        bam_dir = path.join(OutDir, 'bam_dir')
        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       mapping_dir,
                       bam_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class run_star(luigi.Task):
    '''
    run mapping
    '''
    sample = luigi.Parameter()

    def requires(self):
        return prepare()

    def run(self):
        each_sample_mapping_dir = path.join(OutDir, 'mapping_dir', self.sample)
        if not path.exists(each_sample_mapping_dir):
            system('mkdir -p {}'.format(each_sample_mapping_dir))

        tmp = run_cmd(['STAR',
                       '--genomeDir',
                       IndexDir,
                       '--readFilesIn',
                       '{0}/{1}_1.clean.fq.gz'.format(CleanDir, self.sample),
                       '{0}/{1}_2.clean.fq.gz'.format(CleanDir, self.sample),
                       '--readFilesCommand zcat',
                       '--outFileNamePrefix',
                       '{0}/mapping_dir/{1}/'.format(OutDir, self.sample),
                       '--runThreadN',
                       '{}'.format(STAR_THREAD),
                       '--outSAMtype BAM SortedByCoordinate',
                       '--outSAMstrandField intronMotif',
                       '--outFilterType BySJout',
                       '--outFilterMultimapNmax 20',
                       '--alignSJoverhangMin 8',
                       '--alignSJDBoverhangMin 1',
                       '--outFilterMismatchNmax 999',
                       '--alignIntronMin 20',
                       '--alignIntronMax 1000000',
                       '--alignMatesGapMax 1000000'])

        with self.output().open('w') as mapping_log:
            mapping_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.mapping.log'.format(OutDir, self.sample))


class get_bam_file(luigi.Task):
    '''
    link star output bam to bam dir
    '''

    sample = luigi.Parameter()

    def requires(self):
        return [run_star(sample=each_sample) for each_sample in sample_list]

    def run(self):
        link_cmd = ['ln',
                    '-s',
                    '{0}/mapping_dir/{1}/Aligned.sortedByCoord.out.bam'.format(
                        OutDir, self.sample),
                    '{0}/bam_dir/{1}.bam'.format(OutDir, self.sample)]

        index_cmd = ['samtools',
                     'index',
                     '{0}/bam_dir/{1}.bam'.format(OutDir, self.sample)]
        cmd_list = [link_cmd, index_cmd]

        tmp = run_cmd(cmd_list)

        with self.output().open('w') as get_bam_file_log:
            get_bam_file_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.get_bam_file.log'.format(OutDir, self.sample))


class star_mapping_summary(luigi.Task):

    def requires(self):
        return [get_bam_file(sample=each_sample) for each_sample in sample_list]

    def run(self):

        summary_stats_cmd = ['python',
                             STAR_MAPPING_STATS,
                             SampleInf,
                             '{}/mapping_dir/'.format(OutDir),
                             '{}/mapping_stats'.format(OutDir)]

        summary_plot_cmd = ['Rscript',
                            STAR_MAPPING_STATS_PLOT,
                            '--sample_inf',
                            SampleInf,
                            '--mapping_stats',
                            '{}/mapping_stats.plot'.format(OutDir),
                            '--out_dir',
                            OutDir]

        star_mapping_summary_log_inf = run_cmd(
            [summary_stats_cmd, summary_plot_cmd])

        with self.output().open('w') as star_mapping_summary_log:
            star_mapping_summary_log.write(star_mapping_summary_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/logs/star_mapping_summary.log'.format(OutDir))


class star_mapping_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()
    IndexDir = luigi.Parameter()

    def requires(self):
        global OutDir, SampleInf, CleanDir, IndexDir, sample_list
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        IndexDir = self.IndexDir
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        return star_mapping_summary()

    def run(self):
        ignore_files = ['.ignore', 'logs', 'mapping_dir',
                        'bam_dir', 'mapping_stats.plot', 'Rplots.pdf']

        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
