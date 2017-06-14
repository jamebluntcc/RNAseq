'''
Usage:
    run_kegg.py --compare <compare_name> --blast_out <ko_blast> --species <species_abbr> [--background <kegg_backgroud_sp>] --diff_dir <diff_directory> --out_dir <output_dir>

Options:
    -h --help                           Show this screen.
    --compare=<compare_name>            analysis compare name
    --blast_out=<ko_blast>              all gene blast result.
    --species=<species_abbr>            kegg species abbr.
    --background=[<kegg_backgroud_sp>]  kegg background abbr, default is analysis species.
    --diff_dir=<diff_directory>         differential analysis directory.
    --out_dir=<output_dir>              kegg output directory.
'''

from os import system
from os import path
from glob import glob
from docopt import docopt
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
import python_tools
from RNAseq_lib import EXTRACT_INF_BY_ID
from RNAseq_lib import PATHVIEW
from RNAseq_lib import PATHVIEW_CK
from RNAseq_lib import txt_to_excel


class KEGG_enrich:
    def __init__(self):
        self.all_blast_out = None
        self.species = None
        self.background = None
        self.diff_dir = None
        self.out_dir = None
        self.compare = None

    def check_KOBAS_out(self, kobas_out):
        kobas_out_info = open(kobas_out, 'r').readlines()
        flag = True
        for eachline in kobas_out_info:
            if not eachline.startswith("#") and len(eachline.strip().split('\t')) == 9:
                flag = True
                break
        else:
            flag = False
        return flag

    def treat_KEGG_table(self, kegg_output):
        kegg_out_dir, kegg_out_name = path.split(kegg_output)
        kegg_tmp_file = path.join(kegg_out_dir, 'tmp.%s' % kegg_out_name)
        system('mv %s %s' % (kegg_output, kegg_tmp_file))
        if self.check_KOBAS_out(kegg_tmp_file):
            kegg_out_info = open(kegg_output, 'w')
            with open(kegg_tmp_file, 'r') as kegg_tmp_file_info:
                count = 0
                for eachline in kegg_tmp_file_info:
                    if len(eachline.strip().split('\t')) == 9:
                        if count == 0 and eachline.startswith("#"):
                            kegg_out_info.write(eachline)
                            count += 1
                        elif not eachline.startswith("#"):
                            kegg_out_info.write(eachline)
            kegg_out_info.close()
        python_tools.circ_call_process('rm %s' % (kegg_tmp_file))

    def generate_kobas(self, each_blast_out, kegg_output):
        cmd = 'run_kobas.py -i {each_blast_out}  -t blastout:tab -s {self.species} -d K -o {kegg_output} -b {self.background}'.format(
            **locals())
        return cmd

    def run_kegg_pathview(self, each_diff_file):
        cmd_list = []
        pathway_log_dir = path.join(self.out_dir, 'kegg_pathway_logs')
        python_tools.circ_mkdir_unix(pathway_log_dir)
        each_compare_out_dir = path.join(self.out_dir, self.compare)
        each_diff_file_name = path.basename(each_diff_file)
        each_out_prefix = each_diff_file_name.split(
            '.edgeR.DE_results')[0]
        if 'UP' not in each_out_prefix:
            each_out_prefix = '{}.ALL'.format(each_out_prefix)
        kegg_output = path.join(
            each_compare_out_dir, '%s.kegg.enrichment.txt' % (each_out_prefix))
        pathway_outdir = path.join(
            each_compare_out_dir, '%s.pathway' % each_out_prefix)
        pathview_check_log_file = path.join(
            pathway_log_dir, '%s.log' % (each_out_prefix))
        pathview_cmd = 'python %s --kegg_table %s --blast_out %s --species %s --diff_out %s --out_dir %s' % (
            PATHVIEW, kegg_output, self.all_blast_out, self.species, each_diff_file, pathway_outdir)
        pathview_check_cmd = 'python %s --kegg_table %s --pathway_dir %s --log_file %s' % (
            PATHVIEW_CK, kegg_output, pathway_outdir, pathview_check_log_file)
        python_tools.circ_mkdir_unix(pathway_outdir)
        python_tools.circ_call_process(pathview_cmd)
        python_tools.circ_call_process(pathview_check_cmd)
        cmd_list.extend([pathview_cmd, pathview_check_cmd])
        return cmd_list

    def run_KEGG_enrich(self):
        cmd_list = []
        blast_out_dir = path.join(self.out_dir, 'blast_out')
        python_tools.circ_mkdir_unix(blast_out_dir)
        each_compare_diff_dir = path.join(self.diff_dir, self.compare)
        diff_gene_list = glob(
            '{}/*.diffgenes.txt'.format(each_compare_diff_dir))
        each_compare_out_dir = path.join(self.out_dir, self.compare)
        python_tools.circ_mkdir_unix(each_compare_out_dir)
        for each_diff_file in diff_gene_list:
            each_diff_file_name = path.basename(each_diff_file)
            each_out_prefix = each_diff_file_name.split(
                '.edgeR.DE_results')[0]
            each_diff_inf_prefix = each_out_prefix
            if 'UP' not in each_out_prefix:
                each_diff_inf_prefix = each_out_prefix.split('.')[0]
            each_diff_inf_file = path.join(
                each_compare_diff_dir, '{}.edgeR.DE_results.txt'.format(each_diff_inf_prefix))
            kegg_output = path.join(
                each_compare_out_dir, '%s.kegg.enrichment.txt' % (each_out_prefix))
            each_blast_out = path.join(
                blast_out_dir, '%s.blasttab' % (each_out_prefix))
            extract_each_blast_cmd = 'python %s --id %s --table %s --output %s' % (
                EXTRACT_INF_BY_ID, each_diff_file, self.all_blast_out, each_blast_out)
            kegg_cmd = self.generate_kobas(each_blast_out, kegg_output)
            python_tools.circ_call_process(extract_each_blast_cmd)
            cmd_list.append(extract_each_blast_cmd)
            if path.exists(each_blast_out):
                python_tools.circ_call_process(kegg_cmd)
                cmd_list.append(kegg_cmd)
                if path.exists(kegg_output):
                    self.treat_KEGG_table(kegg_output)
                    txt_to_excel(kegg_output)
                    pathway_cmd = self.run_kegg_pathview(each_diff_inf_file)
                    cmd_list.extend(pathway_cmd)
                else:
                    cmd_list.append(
                        "## {} not exists!".format(kegg_output))
            else:
                cmd_list.append("## {} not exists!".format(each_blast_out))
        return cmd_list


if __name__ == '__main__':
    arguments = docopt(__doc__, version='run kegg v1.0')
    my_kegg_enrich = KEGG_enrich()
    my_kegg_enrich.all_blast_out = arguments['--blast_out']
    my_kegg_enrich.species = arguments['--species']
    if arguments['--background']:
        my_kegg_enrich.background = arguments['--background']
    else:
        my_kegg_enrich.background = my_kegg_enrich.species
    my_kegg_enrich.diff_dir = arguments['--diff_dir']
    my_kegg_enrich.out_dir = arguments['--out_dir']
    my_kegg_enrich.compare = arguments['--compare']
    my_kegg_enrich.run_KEGG_enrich()
