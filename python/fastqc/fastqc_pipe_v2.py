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

    def run(self):
        log_dir = os.path.join(OutDir, 'logs')
        fastqc_outdir = os.path.join(OutDir, 'fastqc_results')
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
                        '{}'.format(OutDir)])
        with self.output().open('w') as qc_summary:
            qc_summary.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/fastqc_general_stats.txt'.format(OutDir))

class gc_plot(luigi.Task):
    '''
    plot gc graph
    '''
    sample = luigi.Parameter()

    def requires(self):
        return fastqc_summary()

    def run(self):
        gc_file = os.path.join(OutDir, 'gc_plot', '{0}.gc.txt'.format(self.sample))
        out_prefix = os.path.join(OutDir, 'gc_plot', self.sample)
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
        sample_list = [each.strip() for each in open(self.SampleInf)]
        return reads_quality_plot()

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
