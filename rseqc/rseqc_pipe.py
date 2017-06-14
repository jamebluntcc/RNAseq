import luigi
from os import path
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
# from RNAseq_lib import No_task
from RNAseq_lib import run_cmd
from RNAseq_lib import READ_DISTRIBUTION_PLOT_PREPARE
from RNAseq_lib import RSEQC_PLOT_R
from RNAseq_lib import rsync_pattern_to_file
from python_tools import write_obj_to_file


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''

    OutDir = luigi.Parameter()

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        read_distribution = path.join(OutDir, 'read_distribution')
        genebody_coverage = path.join(OutDir, 'genebody_coverage')
        inner_distance = path.join(OutDir, 'inner_distance')
        junction_saturation = path.join(OutDir, 'junction_saturation')
        read_duplication = path.join(OutDir, 'read_duplication')
        infer_experiment = path.join(OutDir, 'infer_experiment')

        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       read_distribution,
                       genebody_coverage,
                       inner_distance,
                       junction_saturation,
                       read_duplication,
                       infer_experiment])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class read_distribution(luigi.Task):
    '''
    reads distribution on differient region of genome
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['read_distribution.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-r',
                       BedFile])
        with self.output().open('w') as read_distribution_inf:
            read_distribution_inf.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/read_distribution/{1}.read_distribution.txt'.format(OutDir, self.sample))


class genebody_coverage(luigi.Task):

    '''
    genebody coverage
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return [read_distribution(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['geneBody_coverage.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-r',
                       BedFile,
                       '-o',
                       '{0}/genebody_coverage/{1}'.format(OutDir, self.sample)])
        with self.output().open('w') as genebody_coverage_log:
            genebody_coverage_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.genebody_coverage.log'.format(OutDir, self.sample))


class inner_distance(luigi.Task):
    '''
    inner distance

    '''
    OutDir = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        return [genebody_coverage(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['inner_distance.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-r',
                       BedFile,
                       '-o',
                       '{0}/inner_distance/{1}'.format(OutDir, self.sample)])
        with self.output().open('w') as inner_distance_log:
            inner_distance_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.inner_distance.log'.format(OutDir, self.sample))


class junction_saturation(luigi.Task):

    '''
    junction saturation
    '''

    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return [inner_distance(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['junction_saturation.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-r',
                       BedFile,
                       '-o',
                       '{0}/junction_saturation/{1}'.format(OutDir, self.sample)])
        with self.output().open('w') as junction_saturation_log:
            junction_saturation_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.junction_saturation.log'.format(OutDir, self.sample))


class read_duplication(luigi.Task):

    '''
    read duplication

    '''
    OutDir = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        return [junction_saturation(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['read_duplication.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-o',
                       '{0}/read_duplication/{1}'.format(OutDir, self.sample)])
        with self.output().open('w') as read_duplication_log:
            read_duplication_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.read_duplication.log'.format(OutDir, self.sample))


class infer_experiment(luigi.Task):
    '''
    infer experiment

    '''
    OutDir = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        return [read_duplication(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['infer_experiment.py',
                       '-i',
                       '{0}/{1}.bam'.format(BamDir, self.sample),
                       '-r',
                       BedFile])
        with self.output().open('w') as infer_experiment_log:
            infer_experiment_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/infer_experiment/{1}.infer_experiment.txt'.format(OutDir, self.sample))


class read_duplication_plot_prepare(luigi.Task):
    '''
    collection read duplication plot data
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return [infer_experiment(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['python',
                       READ_DISTRIBUTION_PLOT_PREPARE,
                       SampleInf,
                       '{0}/read_distribution/'.format(OutDir)])

        with self.output().open('w') as read_distribution_plot_prepare_logs:
            read_distribution_plot_prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/read_distribution_plot_prepare.log'.format(OutDir))


class rseqc_plot(luigi.Task):

    '''
    rseqc data plot
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return read_duplication_plot_prepare(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['Rscript',
                       RSEQC_PLOT_R,
                       '--sample_inf',
                       SampleInf,
                       '--read_distribution_dir',
                       '{}/read_distribution'.format(OutDir),
                       '--genebody_cov_dir',
                       '{}/genebody_coverage'.format(OutDir),
                       '--inner_distance_dir',
                       '{}/inner_distance'.format(OutDir),
                       '--reads_duplication_dir',
                       '{}/read_duplication'.format(OutDir)])

        with self.output().open('w') as rseqc_plot_logs:
            rseqc_plot_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/rseqc_plot.log'.format(OutDir))


class rseqc_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    BamDir = luigi.Parameter()
    BedFile = luigi.Parameter()
    # UpstreamTask = luigi.Parameter(default=No_task)

    def requires(self):
        global OutDir, SampleInf, BamDir, BedFile, sample_list, UpstreamTask
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        BamDir = self.BamDir
        BedFile = self.BedFile
        # UpstreamTask = self.UpstreamTask
        sample_list = [each.strip().split()[1]
                       for each in open(self.SampleInf)]
        return rseqc_plot(OutDir=OutDir)

    def run(self):
        ignore_files = ['.ignore', 'logs', 'read_duplication/*.DupRate_plot.*',
                        'read_distribution/read_distribution.summary.txt', 'junction_saturation',
                        'inner_distance/*inner_distance_plot*', 'inner_distance/*inner_distance.txt',
                        'infer_experiment', 'genebody_coverage/*geneBodyCoverage.curves.pdf',
                        'genebody_coverage/*geneBodyCoverage.r', 'Rplots.pdf']
        pdf_report_files_pattern = ['inner_distance/*inner_distance.bar.png', 'read_duplication/*reads_duplication.point.png',
                            'genebody_coverage/*genebody_coverage.point.png', 'read_distribution/read_distribution.bar.png', 'read_distribution/*read_distribution.pie.png']
        pdf_report_files = rsync_pattern_to_file(self.OutDir, pdf_report_files_pattern)
        pdf_report_ini = path.join(self.OutDir, '.report_files')
        write_obj_to_file(pdf_report_files, pdf_report_ini)
        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
