import subprocess
import luigi
from os import path
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
    def run(self):
        log_dir = path.join(OutDir, 'logs')
        kallisto_dir = path.join(OutDir, 'kallisto')

        tmp = run_cmd(['mkdir',
                        '-p',
                        log_dir,
                        kallisto_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))

class run_kallisto(luigi.Task):
    '''
    run kallisto
    '''

    sample = luigi.Parameter()

    def requires(self):
        return prepare()

    def run(self)
