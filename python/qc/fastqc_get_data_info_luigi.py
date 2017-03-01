import sys
import os
import pandas as pd

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' sample_list qc_dir > out.summary'
    sys.exit(0)

sample_list_file = sys.argv[1]
qc_dir = sys.argv[2]

sample_list = [each.strip() for each in open(sample_list_file)]
sample_info_dict = {}

GC_HEADER = ['#Base', 'A', 'T', 'G', 'C', 'N-Count']

def add_gc_and_n_inf(gc_info_dict,fastqc_line, fq_number, omit_base = False):
    if fastqc_line.startswith('#'):
        global gc_col_name_list
        gc_col_name_list = fastqc_line.strip().split('\t')
    else:
        eachline_inf = fastqc_line.strip().split('\t')
        for m, each_item in enumerate(eachline_inf):
            each_col_name = gc_col_name_list[m]
            if m == 0:
                each_item = int(each_item.split('-')[0])+(fq_number-1)*sample_info_dict[each_sample][0]
                if omit_base:
                    continue
            if each_col_name not in gc_info_dict:
                gc_info_dict[each_col_name] = [each_item]
            else:
                gc_info_dict[each_col_name].append(each_item)

for each_sample in sample_list:
    each_sample_gc_info_dict = {}
    sample_info_dict[each_sample] = [0, 0, 0, []]
    for n in (1, 2):
        each_qc_dir = os.path.join(qc_dir, '%s_%s_clean.fq_fastqc' % (each_sample, n))
        each_qc_file = os.path.join(each_qc_dir, 'fastqc_data.txt')
        with open(each_qc_file) as each_qc_file_info:
            q30_flag = 0
            gc_flag = 0
            n_flag = 0
            for eachline in each_qc_file_info:
                eachline_info = eachline.strip().split('\t')
                if 'Sequence length' in eachline:
                    sample_info_dict[each_sample][0] = int(eachline_info[1])
                if 'Total Sequences' in eachline:
                    sample_info_dict[each_sample][1] += int(eachline_info[1])
                if eachline.startswith("%GC"):
                    gc_content = int(eachline_info[1])
                    sample_info_dict[each_sample][3].append(gc_content)
                if 'Per sequence quality scores' in eachline:
                    q30_flag = 1
                    continue
                if eachline.startswith('>>'):
                    q30_flag = 0
                    gc_flag = 0
                    n_flag = 0
                if not eachline.startswith('#') and q30_flag == 1 and int(float(eachline_info[0])) > 30 :
                    sample_info_dict[each_sample][2] += int(float(eachline_info[1]))
                if 'Per base sequence content' in eachline:
                    gc_flag = 1
                    continue
                if 'Per base N content' in eachline:
                    n_flag = 1
                    continue
                if gc_flag == 1:
                    add_gc_and_n_inf(each_sample_gc_info_dict, eachline, n)
                if n_flag == 1:
                    add_gc_and_n_inf(each_sample_gc_info_dict, eachline, n, True)
    each_sample_gc_info_df = pd.DataFrame(each_sample_gc_info_dict)
    each_sample_gc_info_table = os.path.join(qc_dir, '{0}.gc.txt'.format(each_sample))
    each_sample_gc_info_df.to_csv(each_sample_gc_info_table, sep = '\t', float_format = '%.3f', columns =  GC_HEADER, index = False)


print 'Sample_ID\tReads_number(M)\tReads_length(bp)\tData_size(G)\tQ30(%)\tGC(%)'
for each_sample in sample_list:
    reads_num = round(sample_info_dict[each_sample][1]/float(10**6), 2)
    read_length = sample_info_dict[each_sample][0]
    data_size = round(sample_info_dict[each_sample][1]*sample_info_dict[each_sample][0]/float(1000**3), 2)
    q30 = round(sample_info_dict[each_sample][2]*100/float(sample_info_dict[each_sample][1]), 2)
    gc_content = sum(sample_info_dict[each_sample][3])/float(len(sample_info_dict[each_sample][3]))
    print '%s\t%s\t%s\t%s\t%s\t%s' % (each_sample, reads_num, read_length, data_size, q30, gc_content)
