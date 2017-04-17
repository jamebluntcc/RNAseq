'''
Usage:
    mRNA_pipe_v1.py [--project_dir=<project_dir>]
                    (--project_name=<project_name>)
                    (--fq_dir=<fq_dir>)
                    (--sample_inf=<sample_inf>)
                    (--transcript=<transcript_fa_file>)
                    (--gene2tr=<gene_transcript_id_map>)
                    (--goseq_anno=<goseq_annotation>)
                    (--topgo_anno <togo_annotation>)
                    (--gene_length=<gene_length_file>)
                    (--kegg_abbr=<species_kegg_abbr>)
                    (--kegg_blast=<kegg_blast_annotation>)

Options:
    --help -h                               show this information
    --project_dir=<project_dir>             project directory
    --project_name=<project_name>           project name
    --fq_dir=<fq_dir>                       fastq directory
    --sample_inf=<sample_inf>               sample information: group<tab>sample
    --transcript=<transcript_fa_file>       transcript fasta file
    --gene2tr=<gene_transcript_id_map>      file map transcript id to gene id
    --goseq_anno=<goseq_annotation>         gene go annotation:(format)gene_id,go_id
    --topgo_anno=<topgo_annotation>         gene go annotation:(format)go_id<tab>gene_id<tab>gene_id ...
    --gene_length=<gene_length_file>        gene length file
    --kegg_abbr=<species_kegg_abbr>         species kegg abbr
    --kegg_blast=<kegg_blast_annotation>    kegg blast out
'''

import ConfigParser
from python_tools import circ_mkdir_unix
from docopt import docopt
import time
from os import system
from os import path
from os import listdir
import subprocess
from collections import OrderedDict

ConfigFile = '/home/public/scripts/RNAseq/python/configure.ini'
conf = ConfigParser.ConfigParser()
conf.read(ConfigFile)
fastqc_pipe = conf.get('server34', 'fastqc_pipe')
rseqc_pipe = conf.get('server34', 'rseqc_pipe')
quant_pipe = conf.get('server34', 'quant_pipe')
enrich_pipe = conf.get('server34', 'enrich_pipe')

PIPE_THREAD = 2

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output

def add_pipe_message(monitor_dir, message = 'pipe_begin'):
    system('touch {0}/{1}_{2}'.format(monitor_dir, message, time.strftime("%Y%m%d_%H%M%S", time.localtime())))
    return 0

def cp_dir_with_ignore(from_dir, to_dir):
    ignore_file = path.join(from_dir, '.ignore')
    system("cp -r {0} {1}".format(from_dir, to_dir))
    with open(ignore_file) as ignore_file_inf:
        for eachline in ignore_file_inf:
            eachline = eachline.strip()
            each_new_ignore_file = path.join(to_dir, eachline)
            system("rm -rf {}".format(each_new_ignore_file))


class Workflow:
    def __init__(self):
        self.work_schedule_dict = {'fastqc':self.run_fastqc, 'quant':self.run_quant, 'enrich':self.run_enrich}
        self.work_schedule_list = ['fastqc', 'quant', 'enrich']
        self.work_status = []
        self.work_dir = []
        self.project_dir = None
        self.project_name = None
        self.monitor_dir = None
        self.fq_dir = None
        self.sample_inf = None
        self.bam_dir = None
        self.bed_file = None
        self.transcript = None
        self.gene2tr = None
        self.goseq_anno = None
        self.topgo_anno = None
        self.gene_length = None
        self.kegg_abbr = None
        self.kegg_blast = None

    def run_fastqc(self):
        add_pipe_message(self.monitor_dir, 'fastqc_start')
        fastqc_dir = path.join(self.project_dir, 'fastqc')

        run_cmd(['python',
        fastqc_pipe,
        'fastqc_collection',
        '--SampleInf',
        self.sample_inf,
        '--OutDir',
        fastqc_dir,
        '--CleanDir',
        self.fq_dir,
        '--workers',
        '{}'.format(PIPE_THREAD)])

        return fastqc_dir

    def run_quant(self):
        add_pipe_message(self.monitor_dir, 'quant_start')
        quant_dir = path.join(self.project_dir, 'quantification')

        run_cmd(['python',
        quant_pipe,
        'quant_collection',
        '--SampleInf',
        self.sample_inf,
        '--CleanDir',
        self.fq_dir,
        '--Transcript',
        self.transcript,
        '--Gene2Tr',
        self.gene2tr,
        '--OutDir',
        quant_dir,
        '--workers',
        '{}'.format(PIPE_THREAD)])

        return quant_dir

    def run_enrich(self):
        add_pipe_message(self.monitor_dir, 'enrich_start')
        enrich_dir = path.join(self.project_dir, 'enrichment')

        run_cmd(['python',
        enrich_pipe,
        'enrichment_collection',
        '--QuantDir',
        '{}/quantification'.format(self.project_dir),
        '--GoseqAnno',
        self.goseq_anno,
        '--TopgoAnno',
        self.topgo_anno,
        '--GeneLen',
        self.gene_length,
        '--KEGGAbbr',
        self.kegg_abbr,
        '--KEGGBlast',
        self.kegg_blast,
        '--OutDir',
        enrich_dir,
        '--workers',
        '{}'.format(PIPE_THREAD)])

        return enrich_dir

    def run_result(self):
        add_pipe_message(self.monitor_dir, 'generate_result_start')
        result_dir = path.join(self.project_dir, self.project_name)
        if path.exists(result_dir):
            system('rm -rf {}'.format(result_dir))
        circ_mkdir_unix(result_dir)
        for each_dir in self.work_dir:
            each_dir_ignore = path.join(each_dir, '.ignore')
            each_dir_name = path.basename(each_dir)
            each_dir_result = path.join(result_dir, each_dir_name)
            cp_dir_with_ignore(each_dir, each_dir_result)

    def run_pipe(self):
        for each_module in self.work_schedule_list:
            if each_module not in self.work_status:
                self.work_dir.append(self.work_schedule_dict[each_module]())
        self.run_result()
        add_pipe_message(self.monitor_dir, 'project_finished')

if __name__ == '__main__':
    ## read arguments
    arguments = docopt(__doc__, version = 'mRNA pipeline v1')
    if not arguments['--project_dir']:
        project_dir = path.abspath('.')
    else:
        project_dir = arguments['--project_dir']
    monitor_dir = path.join(project_dir, 'monitor')
    my_workflow = Workflow()
    my_workflow.project_dir = project_dir
    my_workflow.project_name = arguments['--project_name']
    my_workflow.monitor_dir = monitor_dir
    my_workflow.sample_inf = arguments['--sample_inf']
    my_workflow.fq_dir = arguments['--fq_dir']
    my_workflow.transcript = arguments['--transcript']
    my_workflow.gene2tr = arguments['--gene2tr']
    my_workflow.goseq_anno = arguments['--goseq_anno']
    my_workflow.topgo_anno = arguments['--topgo_anno']
    my_workflow.gene_length = arguments['--gene_length']
    my_workflow.kegg_abbr = arguments['--kegg_abbr']
    my_workflow.kegg_blast = arguments['--kegg_blast']
    map(circ_mkdir_unix, [project_dir, monitor_dir])

    ## running pipeline
    my_workflow.run_pipe()
