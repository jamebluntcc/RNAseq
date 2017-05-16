import subprocess
import luigi
from os import path
from os import system
import sys
from kobas import config

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from RNAseq_lib import PICARD_PATH
from RNAseq_lib import BIOMART_DOWNLOAD
from RNAseq_lib import TOPGO_FORMAT
from RNAseq_lib import TRANSCRIPT_FEATURE
from RNAseq_lib import GO_ANNO
from RNAseq_lib import KEGG_BLAST_TR_TO_GENE
from RNAseq_lib import KEGG_ANNO_EXTRACT
from RNAseq_lib import KEGG_ORGANISM_JSON
from RNAseq_lib import get_kegg_biomart_id
from python_tools import circ_mkdir_unix

STAR_THREAD = '8'
BLAST_THREAD = '8'

## TODO
class down_load_file(luigi.Task):

    def requires(self):
        pass

    def run(self):
        pass

    def output(self):
        pass


class fa_index(luigi.Task):

    def run(self):
        samtools_index_cmd = ['samtools',
        'faidx',
        ref_fa]

        ref_fa_name = path.basename(ref_fa)
        ref_fa_prefix = path.splitext(ref_fa_name)[0]

        picard_index_cmd = ['java',
        '-jar',
        PICARD_PATH,
        'CreateSequenceDictionary',
        'R={}'.format(ref_fa),
        'O={0}/{1}.dict'.format(genome_dir, ref_fa_prefix)]

        cmd_list = [samtools_index_cmd, picard_index_cmd]
        tmp = run_cmd(cmd_list)

        with self.output().open('w') as fa_index_log:
            fa_index_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/fa_index.log'.format(log_dir))

class star_index(luigi.Task):

    def run(self):
        star_index_dir = path.join(annotation_dir, 'star_index')
        circ_mkdir_unix(star_index_dir)

        star_index_cmd = ['STAR',
        '--runThreadN',
        STAR_THREAD,
        '--runMode',
        'genomeGenerate',
        '--sjdbOverhang',
        '149',
        '--genomeFastaFiles',
        ref_fa,
        '--sjdbGTFfile',
        ref_gtf,
        '--genomeDir',
        star_index_dir]

        star_index_log_inf = run_cmd(star_index_cmd)

        with self.output().open('w') as star_index_log:
            star_index_log.write(star_index_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/star_index.log'.format(log_dir))


class transcript_inf_table(luigi.Task):

    def run(self):
        transcript_feature_cmd = ['python',
        TRANSCRIPT_FEATURE,
        '--gtf',
        ref_gtf,
        '--species',
        species_latin,
        '--out_dir',
        annotation_dir]

        transcript_inf_table_log_inf = run_cmd(transcript_feature_cmd)

        with self.output().open('w') as transcript_inf_table_log:
            transcript_inf_table_log.write(transcript_inf_table_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/transcript_inf.log'.format(log_dir))


class go_annotation(luigi.Task):

    def requires(self):
        return transcript_inf_table()

    def run(self):

        download_biomart_go_cmd = ['Rscript',
        BIOMART_DOWNLOAD,
        '--gene_tr_file',
        '{0}/{1}.gene_trans_map.txt'.format(annotation_dir, species_latin),
        '--output',
        '{0}/{1}'.format(annotation_dir, species_latin),
        '--species',
        '{}'.format(species_ensembl)]

        get_topgo_go_cmd = ['python',
        TOPGO_FORMAT,
        '--biomart_go',
        '{0}/{1}.go.txt'.format(annotation_dir, species_latin),
        '--out_dir',
        '{}'.format(annotation_dir)]

        gene_go_anno_cmd = ['python',
        GO_ANNO,
        '{0}/{1}.go.txt'.format(annotation_dir, species_latin),
        '{0}/{1}.go_detail.txt'.format(annotation_dir, species_latin),
        '{0}/{1}.go_anno.txt'.format(annotation_dir, species_latin)]

        go_cmd_list = [download_biomart_go_cmd, get_topgo_go_cmd, gene_go_anno_cmd]
        go_annotation_log_inf = run_cmd(go_cmd_list)

        with self.output().open('w') as go_annotation_log:
            go_annotation_log.write(go_annotation_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/go_annotation.log'.format(log_dir))


class tr_blast(luigi.Task):

    def run(self):

        get_tr_fasta = ['gffread',
        ref_gtf,
        '-g',
        ref_fa,
        '-w',
        '{0}/{1}.transcript.fa'.format(annotation_dir, species_latin)]

        blast_2_db = ['blastx',
        '-query',
        '{0}/{1}.transcript.fa'.format(annotation_dir, species_latin),
        '-db',
        '{0}/{1}.pep.fasta'.format(ko_pep_dir, species_kegg),
        '-evalue',
        '1e-5',
        '-outfmt',
        '6',
        '-max_target_seqs',
        '1',
        '-num_threads',
        BLAST_THREAD,
        '-out',
        '{0}/{1}.tr.kegg.blasttab'.format(annotation_dir, species_latin)]

        tr_blast_cmd_list = [get_tr_fasta, blast_2_db]
        tr_blast_log_inf = run_cmd(tr_blast_cmd_list)

        with self.output().open('w') as tr_blast_log:
            tr_blast_log.write(tr_blast_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/tr_blast.log'.format(log_dir))


class ko_annotation(luigi.Task):

    def requires(self):
        return tr_blast()

    def run(self):

        blast_tr2gene = ['python',
        KEGG_BLAST_TR_TO_GENE,
        ref_gtf,
        '{0}/{1}.tr.kegg.blasttab'.format(annotation_dir, species_latin),
        '{0}/{1}.gene.kegg.blasttab'.format(annotation_dir, species_latin)]

        kobas_anno = ['annotate.py',
        '-i',
        '{0}/{1}.gene.kegg.blasttab'.format(annotation_dir, species_latin),
        '-t',
        'blastout:tab',
        '-s',
        species_kegg,
        '-o',
        '{0}/{1}.gene.ko.anno'.format(annotation_dir, species_latin)]

        extract_ko_anno = ['python',
        KEGG_ANNO_EXTRACT,
        '{0}/{1}.gene.ko.anno'.format(annotation_dir, species_latin),
        '{0}/{1}.gene.ko.anno'.format(annotation_dir, species_latin)]

        kegg_cmd_list = [get_tr_fasta, blast_2_db, blast_tr2gene, kobas_anno, extract_ko_anno]
        kegg_annotation_log_inf = run_cmd(kegg_cmd_list)

        with self.output().open('w') as kegg_annotation_log:
            kegg_annotation_log.write(kegg_annotation_log_inf)

    def output(self):
        return luigi.LocalTarget('{}/kegg_annotation.log'.format(log_dir))


class prepare_collection(luigi.Task):

    ref_fa = luigi.Parameter()
    ref_gtf = luigi.Parameter()
    species_latin = luigi.Parameter()

    def requires(self):
        global ref_fa, genome_dir, ref_gtf, annotation_dir, species_latin, species_ensembl, species_kegg, log_dir, ko_pep_dir
        ref_fa, ref_gtf, species_latin= self.ref_fa, self.ref_gtf, self.species_latin
        species_kegg, species_ensembl = get_kegg_biomart_id(KEGG_ORGANISM_JSON, species_latin)
        genome_dir = path.dirname(ref_fa)
        annotation_dir = path.dirname(ref_gtf)
        log_dir = path.join(annotation_dir, 'logs')
        kobasrc = config.getrc()
        ko_pep_dir = kobasrc['blastdb']
        circ_mkdir_unix(log_dir)
        print ko_pep_dir

        return [fa_index(), star_index(), go_annotation(), ko_annotation()]

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
