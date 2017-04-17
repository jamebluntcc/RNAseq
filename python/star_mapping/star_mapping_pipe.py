import subprocess
import luigi
import os
import sys
import ConfigParser

ConfigFile = '/media/zxchen/3dddd6c4-2700-41a5-a677-15b165fa4e64/scripts/python/configure.ini'
conf = ConfigParser.ConfigParser()
conf.read(ConfigFile)


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
        mapping_dir = os.path.join(self.OutDir, 'mapping_dir')
        bam_dir = os.path.join(self.OutDir, 'bam_dir')
        mapping_summary_dir = os.path.join(self.OutDir, 'mapping_stats')
        tmp = run_cmd(['mkdir',
                        '-p',
                        log_dir,
                        mapping_dir,
                        bam_dir,
                        mapping_summary_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

    def output(self):
        return luigi.LocalTarget('{log_dir}/logs/prepare.log'.format(**locals()))

class star_mapping_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()
    IndexDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir = self.OutDir)

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
