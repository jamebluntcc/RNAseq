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

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output


## run fastqc
class run_fastqc(luigi.Task):
    '''
    run fastqc
    '''
    sample = luigi.Parameter()
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/{self.sample}.log'.format(**locals()))

    def run(self):
        tmp = run_cmd([FASTQC,
                        '{self.CleanDir}/{self.sample}_1.clean.fq.gz'.format(**locals()),
                        '{self.CleanDir}/{self.sample}_2.clean.fq.gz'.format(**locals()),
                        '-o',
                        '{self.OutDir}'.format(**locals())])
        with self.output().open('w') as qc_logs:
            qc_logs.write(tmp)

class qc_report(luigi.Task) :
    '''
    generate qc report

    '''
    ## read configure file
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    #SampleInf = '{}'.format(SampleInf)
    #ConfigFile = luigi.Parameter()
    #ConfigFile = '{}'.format(ConfigFile)

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
        return luigi.LocalTarget('{self.OutDir}/qc.report.txt'.format(**locals()))

if __name__ == '__main__':
    luigi.run()
