import luigi
import socket
import sys
import subprocess
from ConfigParser import ConfigParser
from os import path
from python_tools import load_fn_to_obj


# class No_task(luigi.Task):
#     '''
#     empty task
#     '''
#
#     def run(self):
#         pass
#
#     def output(self):
#         pass

def get_kegg_biomart_id(kegg_map_json, sp_latin):
    kegg_map_dict = load_fn_to_obj(kegg_map_json)
    kegg_sp, biomart_sp = '', ''
    if sp_latin in kegg_map_dict:
        kegg_sp = kegg_map_dict[sp_latin]
    sp_latin_list = sp_latin.split('_')
    biomart_sp = '{0}{1}_gene_ensembl'.format(sp_latin_list[0][0], sp_latin_list[1])
    return kegg_sp, biomart_sp

def get_ip_address():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]

def run_cmd(cmd_obj):
    if not cmd_obj:
        sys.exit('empty cmd')
    elif isinstance(cmd_obj[0], str):
        p = subprocess.Popen(cmd_obj, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
        ret_code = p.wait()
        output = p.communicate()[0]
    elif isinstance(cmd_obj[0], list):
        output_list = []
        for each_cmd in cmd_obj:
            p = subprocess.Popen(each_cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
            ret_code = p.wait()
            output_list.append(p.communicate()[0])
        output = '\n'.join(output_list)
    else:
        sys.exit('unknown cmd format')
    return output

## read confiture
script_path = path.dirname(path.abspath(__file__))
configure_file = path.join(script_path, 'configure.ini')
conf = ConfigParser()
conf.read(configure_file)

## get server information
ip_addr = get_ip_address()
server_name = conf.get('host', ip_addr)

###############
## SOFTWARES ##
###############

## F
FASTQC = conf.get(server_name, 'fastqc')

## G
GATK_PATH = conf.get(server_name, 'gatk_path')

## P
PICARD_PATH = conf.get(server_name, 'picard_path')


#############
## SCRIPTS ##
#############

## B
BIOMART_DOWNLOAD = conf.get(server_name, 'biomart_download')

## E
ENRICH_BARPLOT_R = conf.get(server_name, 'enrich_barplot_r')
EXTRACT_INF_BY_ID = conf.get(server_name, 'extract_inf_by_id')

##  F
FASTQC_PIPE = conf.get(server_name, 'fastqc_pipe')
FASTQC_SUMMERY = conf.get(server_name, 'fastqc_data_info')

## G
GC_PLOT = conf.get(server_name, 'gc_plot')
GO_ANALYSIS_R = conf.get(server_name, 'go_analysis_r')
GO_ANNO = conf.get(server_name, 'go_anno')

## K
KALLISTO_TO_DIFF = conf.get(server_name, 'kallisto_to_diff')
KEGG_ANALYSIS_PYTHON = conf.get(server_name, 'kegg_analysis_python')
KEGG_ANNO_EXTRACT = conf.get(server_name, 'kegg_analysis_python')
KEGG_BLAST_TR_TO_GENE = conf.get(server_name, 'kegg_analysis_python')

## P
PATHVIEW = conf.get(server_name, 'pathview')
PATHVIEW_CK = conf.get(server_name, 'pathview_ck')

## R
READ_DISTRIBUTION_PLOT_PREPARE = conf.get(server_name, 'read_distribution_plot_prepare')
RQ_PLOT = conf.get(server_name, 'rq_plot')
RSEQC_PLOT_R = conf.get(server_name, 'rseqc_plot_r')


## S
STAR_MAPPING_STATS = conf.get(server_name, 'star_mapping_stats')
STAR_MAPPING_STATS_PLOT = conf.get(server_name, 'star_mapping_stats_plot')

## T
TOPGO_FORMAT = conf.get(server_name, 'topgo_format')
TRANSCRIPT_FEATURE = conf.get(server_name, 'transcript_feature')

##########
## DATA ##
##########

## K
KEGG_ORGANISM_TXT = conf.get(server_name, 'kegg_organism_txt')
KEGG_ORGANISM_JSON = conf.get(server_name, 'kegg_organism_json')

## S
SWISSPROT_FASTA = conf.get(server_name, 'swissprot_fasta')
