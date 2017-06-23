from __future__ import division
import socket
import sys
import subprocess
from ConfigParser import ConfigParser
from os import path
from os import system
from os import listdir
from python_tools import load_fn_to_obj
from python_tools import circ_mkdir_unix
import pandas as pd
from glob import glob
from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, Sequence, DateTime
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import and_
from PIL import Image
import datetime


def get_ip_address():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]


# read confiture
script_path = path.dirname(path.abspath(__file__))
configure_file = path.join(script_path, 'configure.ini')
conf = ConfigParser()
conf.read(configure_file)

# get server information
ip_addr = get_ip_address()
server_name = conf.get('host', ip_addr)

###############
#  SOFTWARES  #
###############

# F
FASTQC = conf.get(server_name, 'fastqc')

# G
GATK_PATH = conf.get(server_name, 'gatk_path')

# P
PICARD_PATH = conf.get(server_name, 'picard_path')


#############
#  SCRIPTS  #
#############

# B
BIOMART_DOWNLOAD = conf.get(server_name, 'biomart_download')

# D
DIFF_ANALYSIS = conf.get(server_name, 'diff_analysis')

# E
ENRICH_BARPLOT_R = conf.get(server_name, 'enrich_barplot_r')
EXTRACT_INF_BY_ID = conf.get(server_name, 'extract_inf_by_id')

# F
FASTQC_PIPE = conf.get(server_name, 'fastqc_pipe')
FASTQC_SUMMERY = conf.get(server_name, 'fastqc_data_info')

# G
GET_AS_SUMMARY_PLOT_DATA = conf.get(server_name, 'get_as_summary_plot_data')
GC_PLOT_R = conf.get(server_name, 'gc_plot_r')
GO_ANALYSIS_R = conf.get(server_name, 'go_analysis_r')
GO_ANNO = conf.get(server_name, 'go_anno')

# K
KALLISTO_TO_DIFF = conf.get(server_name, 'kallisto_to_diff')
KALLISTO_TO_TABLE = conf.get(server_name, 'kallisto_to_table')
KEGG_ANALYSIS_PYTHON = conf.get(server_name, 'kegg_analysis_python')
KEGG_ANNO_EXTRACT = conf.get(server_name, 'kegg_anno_extract')
KEGG_BLAST_TR_TO_GENE = conf.get(server_name, 'kegg_blast_tr_to_gene')

# P
PATHVIEW = conf.get(server_name, 'pathview')
PATHVIEW_CK = conf.get(server_name, 'pathview_ck')

# Q
QUANT_REPORT = conf.get(server_name, 'quant_report')

# R
READ_DISTRIBUTION_PLOT_PREPARE = conf.get(
    server_name, 'read_distribution_plot_prepare')
RQ_PLOT = conf.get(server_name, 'rq_plot')
RSEQC_PLOT_R = conf.get(server_name, 'rseqc_plot_r')


# S
SIG_AS_PLOT = conf.get(server_name, 'sig_as_plot')
SNP_PLOT = conf.get(server_name, 'snp_plot')
STAR_MAPPING_STATS = conf.get(server_name, 'star_mapping_stats')
STAR_MAPPING_STATS_PLOT = conf.get(server_name, 'star_mapping_stats_plot')

# T
TOPGO_FORMAT = conf.get(server_name, 'topgo_format')
TRANSCRIPT_FEATURE = conf.get(server_name, 'transcript_feature')

##########
#  DATA  #
##########

# D
DATABASE_DIR = conf.get(server_name, 'database_dir')

# M
MYSQL_DATABASE = conf.get(server_name, 'mysql_database')

# K
KEGG_ORGANISM_TXT = conf.get(server_name, 'kegg_organism_txt')
KEGG_ORGANISM_JSON = conf.get(server_name, 'kegg_organism_json')

# P
PROJECT_DIR = conf.get(server_name, 'project_dir')

# S
SWISSPROT_FASTA = conf.get(server_name, 'swissprot_fasta')

############
#  CONFIG  #
############

REPORT_CFG = conf.get(server_name, 'report_cfg')

Base = declarative_base()


def get_kegg_biomart_id(sp_latin):
    kegg_map_dict = load_fn_to_obj(KEGG_ORGANISM_JSON)
    kegg_sp, biomart_sp = '', ''
    if sp_latin in kegg_map_dict:
        kegg_sp = kegg_map_dict[sp_latin]
    sp_latin_list = sp_latin.split('_')
    biomart_sp = '{0}{1}'.format(sp_latin_list[0][0], sp_latin_list[1])
    return kegg_sp, biomart_sp


class species_annotation_info(Base):
    __tablename__ = 'species_annotation_info'

    id = Column(Integer, Sequence('user_id_seq'),
                primary_key=True, autoincrement=True)

    species_latin = Column(String(100))
    species_database = Column(String(50))
    species_database_version = Column(String(50))
    upload_date = Column(DateTime, default=datetime.datetime.utcnow)

    def __repr__(self):
        return "{0}|{1}|{2}|{3}".format(self.species_latin, self.species_database, self.species_database_version, self.upload_date)


class sepcies_annotation_path:
    def __init__(self):
        self.sp_latin = None
        self.sp_database = 'ensembl'
        self.sp_db_version = None
        self.kegg_abbr = None
        self.genome_fa = None
        self.gtf = None
        self.bedfile = None
        self.star_index = None
        self.transcript = None
        self.gene_tr = None
        self.goseq_ano = None
        self.topgo_ano = None
        self.gene_len = None
        self.kegg_blast = None

    def get_database_dir(self):
        eng = create_engine(MYSQL_DATABASE)
        Session = sessionmaker(bind=eng)
        session = Session()
        if self.sp_db_version:
            rs = session.query(species_annotation_info).filter(and_(species_annotation_info.species_latin == self.sp_latin,
                                                                    species_annotation_info.species_database == self.sp_database,
                                                                    species_annotation_info.species_database_version == self.sp_db_version)).order_by(species_annotation_info.upload_date).all()[-1]
        else:
            rs = session.query(species_annotation_info).filter(and_(species_annotation_info.species_latin == self.sp_latin,
                                                                    species_annotation_info.species_database == self.sp_database)).order_by(species_annotation_info.upload_date).all()[-1]
            self.sp_db_version = str(rs).split('|')[-2]
        sp_database_dirs = glob(
            '{0}/{1}/*/{2}/annotation/{3}'.format(DATABASE_DIR, self.sp_database, self.sp_latin, self.sp_db_version))
        if not sp_database_dirs:
            sys.exit('database [{0}|{1}|{2}] not prepared!'.format(
                self.sp_latin, self.sp_database, self.sp_db_version))
        sp_database_dir = sp_database_dirs[0]
        return sp_database_dir

    def get_anno_inf(self, sp_database_dir=''):
        if not sp_database_dir:
            sp_database_dir = self.get_database_dir()

        def get_annotation_path(x): return path.join(
            sp_database_dir, '{0}.{1}'.format(self.sp_latin, x))
        self.kegg_abbr = get_kegg_biomart_id(self.sp_latin)[0]
        if not self.kegg_abbr:
            self.kegg_abbr = 'ko'
        self.genome_fa = get_annotation_path('genome.fa')
        self.geneme_fai = get_annotation_path('genome.fa.fai')
        self.gtf = get_annotation_path('genome.gtf')
        self.bedfile = get_annotation_path('genome.bed')
        self.star_index = path.join(sp_database_dir, 'star_index')
        self.transcript = get_annotation_path('transcript.fa')
        self.gene_tr = get_annotation_path('gene_trans_map.txt')
        self.goseq_ano = get_annotation_path('go.txt')
        self.topgo_ano = get_annotation_path('go_gene_go.txt')
        self.gene_len = get_annotation_path('gene_length.txt')
        self.kegg_blast = get_annotation_path('gene.kegg.blasttab')


def run_cmd(cmd_obj, is_shell_cmd=False):
    if not cmd_obj:
        sys.exit('empty cmd')
    elif isinstance(cmd_obj[0], str):
        p = subprocess.Popen(
            cmd_obj, universal_newlines=True, stdout=subprocess.PIPE, shell=is_shell_cmd)
        output = p.communicate()[0]
    elif isinstance(cmd_obj[0], list):
        output_list = []
        for each_cmd in cmd_obj:
            # for test
            # try:
            #     p = subprocess.Popen(each_cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
            #     ret_code = p.wait()
            #     output_list.append(p.communicate()[0])
            # except:
            #     print each_cmd
            # for run
            p = subprocess.Popen(
                each_cmd, universal_newlines=True, stdout=subprocess.PIPE, shell=is_shell_cmd)
            output_list.append(p.communicate()[0])
        output = '\n'.join(output_list)
    else:
        sys.exit('unknown cmd format')
    return output


def get_diff_splicing_table(rmats_output, out_dir, pvalue=0.05):
    circ_mkdir_unix(out_dir)
    rmats_output_name = path.basename(rmats_output)
    rmats_output_treat = path.join(out_dir, rmats_output_name)
    system(
        "sed -re 's/\"//g' {0} > {1}".format(rmats_output, rmats_output_treat))
    diff_rmats_out = path.join(out_dir, 'diff.{}'.format(rmats_output_name))
    rmats_output_df = pd.read_table(rmats_output_treat, sep='\t')
    diff_rmats_output_df = rmats_output_df[rmats_output_df.FDR <= 0.05]
    diff_rmats_output_df.to_csv(
        diff_rmats_out, sep='\t', index=False, na_rep='NA')


def get_diff_as_plot_cmd(rmats_results, compare_list, bam_file_list, as_type, out_dir, group_cf):

    plot_cmd = ['rmats2sashimiplot',
                '--b1',
                bam_file_list[0],
                '--b2',
                bam_file_list[1],
                '-t',
                as_type,
                '-e',
                rmats_results,
                '--l1',
                compare_list[0],
                '--l2',
                compare_list[1],
                '-o',
                out_dir,
                '--group-info',
                group_cf]

    return plot_cmd


def rsync_pattern_to_file(from_dir, pattern_list):
    pattern_path_list = [
        '{0}/{1}'.format(from_dir, each_pattern) for each_pattern in pattern_list]
    file_path_list = []
    for each_path in pattern_path_list:
        file_path_list.extend(glob(each_path))
    return [each.split('{}/'.format(from_dir))[1] for each in file_path_list]


def get_enrichment_data(enrichment_dir, plots_num=10):
    go_dir = path.join(enrichment_dir, 'go')
    compare_list = listdir(go_dir)
    pathway_plots = []
    for each_compare in compare_list:
        pathway_plots.extend(rsync_pattern_to_file(enrichment_dir, [
                             'kegg/{}/*ALL.pathway/*pathview.png'.format(each_compare)])[:10])
    go_enrich_table = glob(
        '{}/go/*/*.ALL.go.enrichment.txt'.format(enrichment_dir))[0]
    go_enrich_plot = glob(
        '{}/go/*/*go.enrichment.barplot.png'.format(enrichment_dir))[0]
    go_enrich_table_report = path.join(enrichment_dir, 'report.go.table.txt')
    go_enrich_plot_report = path.join(
        enrichment_dir, 'go.enrichment.barplot.png')
    dag_type = ['CC', 'MF', 'BP']
    for each_type in dag_type:
        each_go_dag_plot = glob(
            '{0}/go/*/DAG/ALL.{1}*png'.format(enrichment_dir, each_type))[0]
        each_go_dag_report_plot = path.join(
            enrichment_dir, '{}.GO.DAG.png'.format(each_type))
        system('cp {0} {1}'.format(each_go_dag_plot, each_go_dag_report_plot))
    system('cp {0} {1}'.format(go_enrich_plot, go_enrich_plot_report))
    system('cut -f1-7 {0} > {1}'.format(go_enrich_table,
                                        go_enrich_table_report))
    kegg_enrich_table = glob(
        '{}/kegg/*/*.ALL.kegg.enrichment.txt'.format(enrichment_dir))[0]
    kegg_enrich_plot = glob(
        '{}/kegg/*/*kegg.enrichment.barplot.png'.format(enrichment_dir))[0]
    kegg_enrich_table_report = path.join(
        enrichment_dir, 'report.kegg.table.txt')
    kegg_enrich_plot_report = path.join(
        enrichment_dir, 'kegg.enrichment.barplot.png')
    kegg_pathway_plot = glob(
        '{}/kegg/*/*ALL.pathway/*pathview.png'.format(enrichment_dir))[0]
    kegg_pathway_report_plot = path.join(enrichment_dir, 'kegg.pathview.png')
    system('cut -f1-7 {0} > {1}'.format(kegg_enrich_table,
                                        kegg_enrich_table_report))
    system('cp {0} {1}'.format(kegg_enrich_plot, kegg_enrich_plot_report))
    system('cp {0} {1}'.format(kegg_pathway_plot, kegg_pathway_report_plot))
    repor_data = ['report.go.table.txt', 'report.kegg.table.txt',
                  'go.enrichment.barplot.png', 'kegg.enrichment.barplot.png',
                  'kegg.pathview.png']
    go_dag_plots = ['{}.GO.DAG.png'.format(each) for each in dag_type]
    repor_data.extend(go_dag_plots)
    repor_data.extend(pathway_plots)
    return repor_data


def txt_to_excel(txt_file, sheet_name='Sheet1'):
    pd.formats.format.header_style = None
    txt_df = pd.read_table(txt_file)
    txt_file_name = path.basename(txt_file)
    txt_file_dir = path.dirname(txt_file)
    txt_file_prefix = path.splitext(txt_file_name)[0]
    excel_file = path.join(txt_file_dir, '{}.xlsx'.format(txt_file_prefix))
    writer = pd.ExcelWriter(excel_file, engine='xlsxwriter', options={
                            'strings_to_urls': False})
    txt_df.to_excel(writer, sheet_name, index=False)
    writer.save()


def check_rseqc_condition(genome_fai, longest_chr_size=500000000):
    genome_fai_df = pd.read_table(genome_fai, header=None)
    if max(genome_fai_df.loc[:, 1]) > longest_chr_size:
        return False
    else:
        return True


def add_prefix_to_filename(file_path, add_inf='pdf'):
    file_name = path.basename(file_path)
    file_dir = path.dirname(file_path)
    file_prefix, file_sufix = path.splitext(file_name)
    return path.join(file_dir, '{0}.{1}{2}'.format(file_prefix, add_inf, file_sufix))


# def resize_plot(ori_plot, resize, out_plot):
#     return ['convert', '-resize', '{0}%'.format(resize), ori_plot, out_plot]

def plot_resize(ori_size, target_size):
    factor_size1 = round(target_size[0] / ori_size[0], 2)
    factor_size2 = round(target_size[1] / ori_size[1], 2)
    factor_size = min(factor_size1, factor_size2)
    return [int(each * factor_size) for each in ori_size]


def resize_report_plot(report_dir):
    config = ConfigParser()
    config.read(REPORT_CFG)
    items = config.items('pdf_path')

    for value in items:
        item_name = value[0]
        item_value = value[1]
        each_plot_path = path.join(report_dir, item_value)
        if path.exists(each_plot_path):
            each_plot_resize_path = add_prefix_to_filename(each_plot_path)
            each_plot_img = Image.open(each_plot_path)
            each_plot_img_size = each_plot_img.size
            each_plot_pdf_size = config.get('pdf_size', item_name)
            each_plot_pdf_size = [int(each)
                                  for each in each_plot_pdf_size.split(',')]
            new_size = plot_resize(each_plot_img_size, each_plot_pdf_size)
            each_plot_img.resize(new_size).save(each_plot_resize_path)

    return 'finished resize plot'


def main():
    pass


if __name__ == '__main__':
    main()
