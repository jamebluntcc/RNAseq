import subprocess
import luigi
import os
import sys
import ConfigParser

ConfigFile = '/media/zxchen/3dddd6c4-2700-41a5-a677-15b165fa4e64/scripts/python/configure.ini'
conf = ConfigParser.ConfigParser()
conf.read(ConfigFile)
FASTQC_SUMMERY = conf.get('server34', 'fastqc_data_info')
FASTQC = conf.get('server34', 'fastqc')
GC_PLOT = conf.get('server34', 'gc_plot')
RQ_PLOT = conf.get('server34', 'rq_plot')

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''
    OutDir = luigi.Parameter()

    def run(self):
        log_dir = os.path.join(self.OutDir, 'logs')
        fastqc_outdir = os.path.join(self.OutDir, 'fastqc_results')
        tmp = run_cmd(['mkdir',
                        '-p',
                        log_dir,
                        fastqc_outdir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

    def output(self):
        log_dir = os.path.join(self.OutDir, 'logs')
        return luigi.LocalTarget('{log_dir}/prepare.log'.format(**locals()))


## run fastqc
class run_fastqc(luigi.Task):
    '''
    run fastqc
    '''
    sample = luigi.Parameter()
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir = self.OutDir)

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/logs/{self.sample}.log'.format(**locals()))

    def run(self):
        tmp = run_cmd([FASTQC,
                        '{self.CleanDir}/{self.sample}_1.clean.fq.gz'.format(**locals()),
                        '{self.CleanDir}/{self.sample}_2.clean.fq.gz'.format(**locals()),
                        '--extract',
                        '-o',
                        '{self.OutDir}/fastqc_results'.format(**locals())])
        with self.output().open('w') as qc_logs:
            qc_logs.write(tmp)

class fastqc_summary(luigi.Task) :
    '''
    generate summary table
    '''
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()

    def requires(self):
        sample_list = [each.strip() for each in open(self.SampleInf)]
        return [run_fastqc(sample = sample, CleanDir = self.CleanDir, OutDir = self.OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['python',
                        FASTQC_SUMMERY,
                        '{}'.format(self.SampleInf),
                        '{self.OutDir}'.format(**locals())])
        with self.output().open('w') as qc_summary:
            qc_summary.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/fastqc_general_stats.txt'.format(**locals()))

class gc_plot(luigi.Task):
    '''
    plot gc graph
    '''

    sample = luigi.Parameter()
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()

    def requires(self):
        return fastqc_summary(CleanDir = self.CleanDir, OutDir = self.OutDir, SampleInf = self.SampleInf)

    def run(self):
        gc_file = os.path.join(self.OutDir, 'gc_plot', '{0}.gc.txt'.format(self.sample))
        out_prefix = os.path.join(self.OutDir, 'gc_plot', self.sample)
        tmp = run_cmd(['Rscript',
        GC_PLOT,
        '--gcfile',
        gc_file,
        '--prefix',
        out_prefix])
        with self.output().open('w') as gc_plot_log:
            gc_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/logs/{self.sample}.gc.log'.format(**locals()))

class reads_quality_plot(luigi.Task):
    '''
    plot reads quality
    '''

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()

    def requires(self):
        sample_list = [each.strip() for each in open(self.SampleInf)]
        return [gc_plot(sample = sample, OutDir = self.OutDir, SampleInf = self.SampleInf, CleanDir = self.CleanDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['Rscript',
        RQ_PLOT,
        '{self.OutDir}/reads_quality_plot/'.format(**locals())])
        with self.output().open('w') as reads_quality_plot_log:
            reads_quality_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/logs/reads_quality_plot.log'.format(**locals()))

class fastqc_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()

    def requires(self):
        return reads_quality_plot(OutDir = self.OutDir, SampleInf = self.SampleInf, CleanDir = self.CleanDir)

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
