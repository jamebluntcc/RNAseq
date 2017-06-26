#! /usr/bin/python

import luigi
from os import path
from RNAseq_lib import run_cmd
from RNAseq_lib import qc_info
import fastqc_pipe_v2 as fastqc_pipe
import star_mapping_pipe_v2 as star_mapping_pipe
import rseqc_pipe as rseqc_pipe


class prepare(luigi.Task):

    '''
    prepare directory and others
    '''

    OutDir = luigi.Parameter()

    def run(self):

        log_dir = path.join(self.OutDir, 'logs')
        qc_data_dir = path.join(self.OutDir, 'qc_data')
        prepare_inf = run_cmd(['mkdir',
                               '-p',
                               log_dir,
                               qc_data_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write(prepare_inf)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(self.OutDir))


class get_qc_data(luigi.Task):
    '''
    cut 5M reads of each sample for qc
    '''

    proj_dir = luigi.Parameter()
    sample = luigi.Parameter()
    clean_dir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=self.proj_dir)

    def run(self):
        qc_data_dir = path.join(self.proj_dir, 'qc_data')
        cut_cmds = []
        fq_files = [
            '{0}/{1}_{2}.clean.fq.gz'.format(self.clean_dir, self.sample, each) for each in (1, 2)]
        qc_fq_files = [
            '{0}/{1}_{2}.clean.fq.gz'.format(qc_data_dir, self.sample, each) for each in (1, 2)]
        for n, each_fq in enumerate(fq_files):
            each_cut_cmd = [
                'zcat {0} | head -20000000 | gzip > {1}'.format(each_fq, qc_fq_files[n])]
            cut_cmds.append(each_cut_cmd)
        cut_cmd_inf = run_cmd(cut_cmds, is_shell_cmd=True)
        with self.output().open('w') as get_qc_data_log:
            get_qc_data_log.write(cut_cmd_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.get_qc_data.log'.format(self.proj_dir, self.sample))


class fastqc(luigi.Task):

    proj_dir = luigi.Parameter()
    fq_dir = luigi.Parameter()
    dir_name = 'fastqc'

    def requires(self):
        return [get_qc_data(sample=each_sample, clean_dir=CleanDir, proj_dir=self.proj_dir) for each_sample in sample_list]

    def run(self):
        out_dir = path.join(self.proj_dir, self.dir_name)
        yield fastqc_pipe.fastqc_collection(OutDir=out_dir, SampleInf=SampleInf, CleanDir=self.fq_dir)
        with self.output().open('w') as fastqc_log_inf:
            fastqc_log_inf.write('fastqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/fastqc.log'.format(log_dir))


class mapping(luigi.Task):

    proj_dir = luigi.Parameter()
    fq_dir = luigi.Parameter()
    dir_name = 'mapping'

    def requires(self):
        return [get_qc_data(sample=each_sample, clean_dir=CleanDir, proj_dir=self.proj_dir) for each_sample in sample_list]

    def run(self):
        out_dir = path.join(self.proj_dir, self.dir_name)
        yield star_mapping_pipe.star_mapping_collection(OutDir=out_dir, IndexDir=StarIndex, SampleInf=SampleInf, CleanDir=self.fq_dir)

        with self.output().open('w') as mapping_log_inf:
            mapping_log_inf.write('mapping finished!')

    def output(self):
        return luigi.LocalTarget('{}/mapping.log'.format(log_dir))


class rseqc(luigi.Task):

    proj_dir = luigi.Parameter()
    dir_name = 'rseqc'

    def requires(self):
        return mapping(proj_dir=self.proj_dir, fq_dir=qc_fq_dir)

    def run(self):
        bam_dir = path.join(self.proj_dir, 'mapping', 'bam_dir')
        out_dir = path.join(self.proj_dir, self.dir_name)
        yield rseqc_pipe.rseqc_collection(OutDir=out_dir, SampleInf=SampleInf, BamDir=bam_dir, BedFile=BedFile)

        with self.output().open('w') as rseqc_log:
            rseqc_log.write('rseqc finished')

    def output(self):
        return luigi.LocalTarget('{}/rseqc.log'.format(log_dir))


class qc_collection(luigi.Task):

    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    StarIndex = luigi.Parameter()
    BedFile = luigi.Parameter()

    def requires(self):
        global SampleInf, CleanDir, OutDir, log_dir, sample_list, BedFile
        global qc_fq_dir, StarIndex
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        OutDir = self.OutDir
        BedFile = self.BedFile
        StarIndex = self.StarIndex
        log_dir = path.join(OutDir, 'logs')
        qc_fq_dir = path.join(OutDir, 'qc_data')
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        return [fastqc(proj_dir=OutDir, fq_dir=qc_fq_dir),
                rseqc(proj_dir=OutDir)]

    def run(self):

        my_qc_info = qc_info(self.SampleInf, self.OutDir)
        my_qc_info.check_data()

        with self.output().open('w') as qc_collection_log:
            qc_collection_log.write('qc finished.')

    def output(self):
        return luigi.LocalTarget('{}/logs/qc_collection.log'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
