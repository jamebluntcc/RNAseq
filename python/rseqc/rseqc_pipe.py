import subprocess
import luigi
import os
import sys
import ConfigParser

ConfigFile = '/home/public/scripts/RNAseq/python/configure.ini'
conf = ConfigParser.ConfigParser()
conf.read(ConfigFile)

read_distribution_plot_prepare = conf.get('server34', 'read_distribution_plot_prepare')
r_home = conf.get('server34', 'r_home')
rseqc_plot_r = conf.get('server34', 'rseqc_plot_r')

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
        read_distribution = os.path.join(OutDir, 'read_distribution')
        genebody_coverage = os.path.join(OutDir, 'genebody_coverage')
        inner_distance = os.path.join(OutDir, 'inner_distance')
        junction_saturation = os.path.join(OutDir, 'junction_saturation')
        read_duplication = os.path.join(OutDir, 'read_duplication')
        infer_experiment = os.path.join(OutDir, 'infer_experiment')

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
            prepare_logs.write('prepare finished')

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class read_distribution(luigi.Task):
    '''
    reads distribution on differient region of genome
    '''

    sample = luigi.Parameter()

    def requires(self):
        return prepare()

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

    def requires(self):
        return [read_distribution(sample = sample) for sample in sample_list]

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

    sample = luigi.Parameter()

    def requires(self):
        return [genebody_coverage(sample = sample) for sample in sample_list]

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

    def requires(self):
        return [inner_distance(sample = sample) for sample in sample_list]

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

    sample = luigi.Parameter()

    def requires(self):
        return [junction_saturation(sample = sample) for sample in sample_list]

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

    sample = luigi.Parameter()

    def requires(self):
        return [read_duplication(sample = sample) for sample in sample_list]

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

    def requires(self):
        return [infer_experiment(sample = sample) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['python',
        read_distribution_plot_prepare,
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

    def requires(self):
        return read_duplication_plot_prepare()

    def run(self):
        tmp = run_cmd(['{}/Rscript'.format(r_home),
        rseqc_plot_r,
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

    def requires(self):
        global OutDir, SampleInf, BamDir, BedFile, sample_list
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        BamDir = self.BamDir
        BedFile = self.BedFile
        sample_list = [each.strip().split()[1] for each in open(self.SampleInf)]
        return rseqc_plot()

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
