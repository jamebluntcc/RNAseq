'''
Usage:
    transcript_feature.py --gtf=<gtf_file> --species=<species_latin_name> --out_dir=<ouput_dir>

Options:
    --help -h                       show this information
    --gtf=<gtf_file>                gtf file
    --species=<species_latin_name>  species latin name
    --out_dir=<ouput_dir>           output directory
'''
from docopt import docopt
from os import path
from os import stat
import sys
from HTSeq import GFF_Reader

script_path = os.path.dirname(os.path.abspath(__file__))
RNAseq_lib_path = os.path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_tools import get_transcript_info
from RNAseq_tools import get_gene_info
from python_tools import Median


def get_target_length_and_pos_table(gtf_file, gene_length, gene_pos):
    target_length_info = open(gene_length,'w')
    pos_file_inf = open(gene_pos, 'w')
    transcript_dict = get_transcript_info(gtf_file)
    gene_dict = get_gene_info(transcript_dict)
    for each_gene in gene_dict :
        each_gene_len = Median(gene_dict[each_gene]['transcript_len'])
        target_length_info.write('{0}\t{1}\n'.format(each_gene, each_gene_len))
        each_start = min(gene_dict[each_gene]['tss'])
        each_end = max(gene_dict[each_gene]['tts'])
        each_chrom = gene_dict[each_gene]['chrom']
        each_strand = gene_dict[each_gene]['strand']
        pos_file_inf.write('{0}\t{1}:{2}-{3}\t{4}\n'.format(each_gene, each_chrom, each_start, each_end, each_strand))
    target_length_info.close()
    pos_file_inf.close()
    return 'produced gene length and position file!'

if __name__ == '__main__':
    ## read arguments
    arguments = docopt(__doc__, version = '1.0')
    gtf_file = arguments['--gtf']
    species = arguments['--species']
    out_dir = arguments['--out_dir']
    ## gene length file and gene locus file
    gene_length_file = path.join(out_dir, '{}.gene_length.txt'.format(species))
    gene_locus_file = path.join(out_dir, '{}.gene_locus.txt'.format(species))
    get_target_length_table(gtf_file, gene_length_file, gene_locus_file)

    ## gene transcript map file
    tr_dict ={}
    gene_tr_map_file = path.join(out_dir,'{}.gene_trans_map.txt'.format(species))
    with open(gene_tr_map_file, 'w') as gene_tr_map_file_inf:
        for eachline in GFF_Reader(gtf_file) :
            if 'transcript_id' in eachline.attr:
                transcript_id = eachline.attr['transcript_id']
                gene_id = eachline.attr['gene_id']
                if transcript_id not in tr_dict:
                    tr_dict[transcript_id] = gene_id
                    gene_tr_map_file_inf.write('{0}\t{1}\n'.format(gene_id,transcript_id))
