import luigi
from os import path
import sys
from glob import glob
from kobas import config

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from python_tools import circ_mkdir_unix
from python_tools import write_obj_to_file
# from RNAseq_lib import SWISSPROT_FASTA

BLAST_THREAD = '2'
STRINGTIE_THREAD = '2'


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        assembly_dir = path.join(OutDir, 'assembly_dir')
        annotation_dir = path.join(OutDir, 'annotation')

        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       assembly_dir,
                       annotation_dir])

        with self.output().open('w') as prepare_logs:
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class stringtie_assembly(luigi.Task):

    sample = luigi.Parameter()

    def requires(self):
        return prepare()

    def run(self):

        stringtie_cmd = ['stringtie',
                         '{0}/{1}.bam'.format(BamDir, self.sample),
                         '-p',
                         STRINGTIE_THREAD,
                         '-o',
                         '{0}/assembly_dir/{1}.gtf'.format(OutDir, self.sample)]

        if RefGtf:
            stringtie_cmd.extend(['-G', RefGtf])

        stringtie_assembly_log_inf = run_cmd(stringtie_cmd)

        with self.output().open('w') as stringtie_assembly_log:
            stringtie_assembly_log.write(stringtie_assembly_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.assembly.log'.format(OutDir, self.sample))


class stringtie_merge(luigi.Task):

    def requires(self):
        return [stringtie_assembly(sample=each_sample) for each_sample in sample_list]

    def run(self):
        gtf_list_file = path.join(OutDir, 'assembly_dir', 'gtf.list')
        gtf_path_list = [path.join(OutDir, 'assembly_dir', '{}.gtf'.format(
            each_sample)) for each_sample in sample_list]
        write_obj_to_file(gtf_path_list, gtf_list_file)

        merge_gtf_cmd = ['stringtie',
                         '--merge',
                         '-m',
                         '200',
                         '-T',
                         '0.1',
                         '-f',
                         '0.1',
                         '-o',
                         '{}/stringtie_merge.gtf'.format(OutDir),
                         gtf_list_file]
        if RefGtf:
            merge_gtf_cmd.extend(['-G', RefGtf])

        get_fa_cmd = ['gffread',
                      '{}/stringtie_merge.gtf'.format(OutDir),
                      '-g',
                      RefFa,
                      '-w',
                      '{}/stringtie_merge.fa'.format(OutDir)]

        stringtie_merge_log_inf = run_cmd(
            [merge_gtf_cmd, get_fa_cmd])

        with self.output().open('w') as stringtie_merge_log:
            stringtie_merge_log.write(stringtie_merge_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/stringtie_merge.log'.format(OutDir))


class blast_annotate(luigi.Task):

    database_name = luigi.Parameter()
    database = luigi.Parameter()
    fasta_file = luigi.Parameter()
    blast_out_dir = luigi.Parameter()

    def requires(self):
        return prepare()

    def run(self):
        fasta_file_name = path.basename(self.fasta_file)

        # orf_pred_cmd = ['TransDecoder.LongOrfs',
        #                 '-t',
        #                 '{}/stringtie_merge.fa'.format(OutDir),
        #                 '--gene_trans_map',
        #                 SampleInf]
        # link_orf_cmd = ['ln',
        #                 '-s',
        #                 '{}/stringtie_merge.fa.transdecoder_dir/longest_orfs.pep'.format(
        #                     OutDir),
        #                 '{}/annotation']

        blast_cmd = ['blastp',
                     '-query',
                     self.fasta_file,
                     '-db',
                     self.database,
                     '-evalue',
                     '1e-5',
                     '-outfmt',
                     '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
                     '-max_target_seqs',
                     '1',
                     '-num_threads',
                     BLAST_THREAD,
                     '-out',
                     '{0}/{1}.blasttab'.format(self.blast_out_dir, fasta_file_name)]

        blast_log_inf = run_cmd(blast_cmd)

        with self.output().open('w') as blast_log:
            blast_log.write(blast_log_inf)

    def output(self):
        fasta_file_name = path.basename(self.fasta_file)
        return luigi.LocalTarget('{0}/logs/{1}_{2}.blast.log'.format(OutDir, self.database_name, fasta_file_name))


class ref_assembly_collection(luigi.Task):

    SampleInf = luigi.Parameter()
    BamDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    RefGtf = luigi.Parameter(default='')
    RefFa = luigi.Parameter()
    # KEGGAbbr = luigi.Parameter(default='')

    def requires(self):
        global SampleInf, BamDir, sample_list, RefGtf, RefFa, OutDir
        SampleInf = self.SampleInf
        BamDir = self.BamDir
        sample_list = [each.strip().split()[1]
                       for each in open(self.SampleInf)]
        RefGtf = self.RefGtf
        RefFa = self.RefFa
        OutDir = path.abspath(self.OutDir)
        # KEGG blast
        # KEGGAbbr = self.KEGGAbbr
        # kobasrc = config.getrc()
        # ko_pep_dir = kobasrc['blastdb']
        # ko_pep = path.join(ko_pep_dir, '{}.pep.fasta'.format(KEGGAbbr))
        # ko_out_dir = path.join(OutDir, 'kegg')
        # circ_mkdir_unix(ko_out_dir)
        # return [blast_annotate(database_name='kegg', database=ko_pep,
        # fasta_file=each_fa, blast_out_dir=ko_out_dir) for each_fa in
        # split_fastas]
        return stringtie_merge()

        # SWISS-PROT blast
        # swissprot_out_dir = path.join(OutDir, 'swissprot')
        # circ_mkdir_unix(swissprot_out_dir)

        # return [blast_annotate(database_name = 'swissprot', database =
        # SWISSPROT_FASTA, fasta_file = each_fa, blast_out_dir =
        # swissprot_out_dir) for each_fa in split_fastas]

    def run(self):
        pass

    def output(self):
        pass


if __name__ == '__main__':
    luigi.run()
