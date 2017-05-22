#! /usr/bin/python

import luigi
from os import path
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import FASTQC_SUMMERY
from RNAseq_lib import FASTQC
from RNAseq_lib import GC_PLOT
from RNAseq_lib import RQ_PLOT
from RNAseq_lib import run_cmd

class prepare(luigi.Task):
    '''
    prepare directory and others
    '''

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        fastqc_outdir = path.join(OutDir, 'fastqc_results')
        tmp = run_cmd(['mkdir',
                        '-p',
                        log_dir,
                        fastqc_outdir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))

## run fastqc
class run_fastqc(luigi.Task):
    '''
    run fastqc
    '''
    sample = luigi.Parameter()

    def requires(self):
        return prepare()

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.log'.format(OutDir, self.sample))

    def run(self):
        tmp = run_cmd([FASTQC,
                        '{0}/{1}_1.clean.fq.gz'.format(CleanDir, self.sample),
                        '{0}/{1}_2.clean.fq.gz'.format(CleanDir, self.sample),
                        '--extract',
                        '-o',
                        '{}/fastqc_results'.format(OutDir)])
        with self.output().open('w') as qc_logs:
            qc_logs.write(tmp)

class fastqc_summary(luigi.Task) :
    '''
    generate summary table
    '''

    def requires(self):
        return [run_fastqc(sample = sample) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['python',
                        FASTQC_SUMMERY,
                        '{}'.format(SampleInf),
                        '{}'.format(OutDir),
                        '{}/fastqc_general_stats'.format(OutDir)])
        with self.output().open('w') as qc_summary:
            qc_summary.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/fastqc_summary.log'.format(OutDir))

class gc_plot(luigi.Task):
    '''
    plot gc graph
    '''
    sample = luigi.Parameter()

    def requires(self):
        return fastqc_summary()

    def run(self):
        gc_file = path.join(OutDir, 'gc_plot', '{0}.gc.txt'.format(self.sample))
        out_prefix = path.join(OutDir, 'gc_plot', self.sample)
        tmp = run_cmd(['Rscript',
        GC_PLOT,
        '--gcfile',
        gc_file,
        '--prefix',
        out_prefix])
        with self.output().open('w') as gc_plot_log:
            gc_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.gc.log'.format(OutDir, self.sample))

class reads_quality_plot(luigi.Task):
    '''
    plot reads quality
    '''

    def requires(self):
        return [gc_plot(sample = sample) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['Rscript',
        RQ_PLOT,
        '{}/reads_quality_plot/'.format(OutDir)])
        with self.output().open('w') as reads_quality_plot_log:
            reads_quality_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/reads_quality_plot.log'.format(OutDir))

class fastqc_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()

    def requires(self):
        global OutDir, SampleInf, CleanDir, sample_list
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        return reads_quality_plot()

    def run(self):
        ignore_files = ['.ignore', 'logs', 'fastqc_results/*zip']
        #include_files = ['fastqc_results/*zip', 'fastqc_results/*html', 'reads_quality_plot', 'gc_plot', 'fastqc_general_stats.*']
        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))

if __name__ == '__main__':
    luigi.run()
