import json
import re
from os import path
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import KEGG_ORGANISM_TXT
from RNAseq_lib import KEGG_ORGANISM_JSON
from python_tools import write_obj_to_json

kegg_name_map_dict = {}
with open(KEGG_ORGANISM_TXT) as kegg_organism_txt_inf:
    for eachline in kegg_organism_txt_inf:
        eachline_inf = eachline.rstrip().split('\t')
        kegg_sp = eachline_inf[1]
        latin_info = eachline_inf[2]
        if '(' in latin_info:
            latin_name = re.match(r'(.*)\(', latin_info).groups()[0].lower().strip()
            latin_name = re.sub(' ', '_', latin_name)
        else:
            latin_name = re.sub(' ', '_', latin_info.lower())
        kegg_name_map_dict[latin_name] = kegg_sp

write_obj_to_json(kegg_name_map_dict, KEGG_ORGANISM_JSON)
