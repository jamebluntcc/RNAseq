import subprocess
import luigi
from os import path
from os import system
import ConfigParser

STAR_THREAD = 8

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
        mapping_dir = path.join(OutDir, 'mapping_dir')
        bam_dir = path.join(OutDir, 'bam_dir')
        mapping_summary_dir = path.join(OutDir, 'mapping_stats')
        tmp = run_cmd(['mkdir',
                        '-p',
                        log_dir,
                        mapping_dir,
                        bam_dir,
                        mapping_summary_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

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
        '--outFilterType BySJout',
        '--outFilterMultimapNmax 20',
        '--alignSJoverhangMin 8',
        '--alignSJDBoverhangMin 1',
        '--outFilterMismatchNmax 999',
        '--alignIntronMin 20',
        '--alignIntronMax 1000000',
        '--alignMatesGapMax 1000000',
        '--outSAMstrandField intronMotif'])

        with self.output().open('w') as mapping_log:
            mapping_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.mapping.log'.format(OutDir, self.sample))


class star_mapping_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()
    IndexDir = luigi.Parameter()
    ConfigFile = luigi.Parameter()

    def requires(self):
        global OutDir, SampleInf, CleanDir, IndexDir, sample_list, ConfigFile, conf
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        IndexDir = self.IndexDir
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        ConfigFile = self.ConfigFile
        conf = ConfigParser.ConfigParser()
        conf.read(ConfigFile)
        return [run_star(sample = each_sample) for each_sample in sample_list]

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
