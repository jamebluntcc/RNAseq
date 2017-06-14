import sys
import os
import pandas as pd
import xlsxwriter

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' sample_list qc_dir out.summary.prefix'
    sys.exit(0)


def circ_mkdir_unix(path):
    cmd = 'mkdir -p %s' % path
    if not os.path.isdir(path):
        os.system(cmd)


def output_to_xlsx(worksheet, info_list, row_num):
    for n, each in enumerate(info_list):
        each = str(each)
        if isinstance(each, str):
            worksheet.write_string(row_num, n, each)
        else:
            worksheet.write_number(row_num, n, each)


def add_gc_and_n_inf(gc_info_dict, fastqc_line, fq_number, omit_base=False):
    if fastqc_line.startswith('#'):
        global gc_col_name_list
        gc_col_name_list = fastqc_line.strip().split('\t')
    else:
        eachline_inf = fastqc_line.strip().split('\t')
        for m, each_item in enumerate(eachline_inf):
            each_col_name = gc_col_name_list[m].split('-')[0]
            if m == 0:
                each_item = int(each_item.split(
                    '-')[0]) + (fq_number - 1) * sample_info_dict[each_sample][0]
                if omit_base:
                    continue
            if each_col_name not in gc_info_dict:
                gc_info_dict[each_col_name] = [each_item]
            else:
                gc_info_dict[each_col_name].append(each_item)


sample_list_file = sys.argv[1]
qc_dir = sys.argv[2]
out_summary_prefix = sys.argv[3]

sample_list = [each.strip().split()[1] for each in open(sample_list_file)]
sample_info_dict = {}

GC_HEADER = ['#Base', 'A', 'T', 'G', 'C', 'N']

reads_quality_dir = os.path.join(qc_dir, 'reads_quality_plot')
gc_dir = os.path.join(qc_dir, 'gc_plot')
circ_mkdir_unix(reads_quality_dir)
circ_mkdir_unix(gc_dir)
# merged_quality_file = os.path.join(reads_quality_dir, 'example.data.txt')
# merged_quality_file_inf = open(merged_quality_file, 'w')
# merged_quality_file_inf.write('Sample_ID\tQuality\tCount\tPercent\n')

for each_sample in sample_list:
    each_sample_gc_info_dict = {}
    sample_info_dict[each_sample] = [0, 0, 0, []]
    reads_quality_dict = {}
    each_reads_quality_file = os.path.join(
        reads_quality_dir, '{}.reads_quality.txt'.format(each_sample))
    for n in (1, 2):
        each_qc_dir = os.path.join(
            qc_dir, 'fastqc_results', '%s_%s.clean_fastqc' % (each_sample, n))
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
                if not eachline.startswith('#') and q30_flag == 1:
                    if int(float(eachline_info[0])) > 30:
                        sample_info_dict[each_sample][2] += int(
                            float(eachline_info[1]))
                    each_quality = int(eachline_info[0])
                    each_count = float(eachline_info[1])
                    if each_sample in reads_quality_dict and each_quality in reads_quality_dict[each_sample]:
                        reads_quality_dict[each_sample][each_quality] += each_count
                    else:
                        reads_quality_dict.setdefault(each_sample, {})[
                            each_quality] = each_count
                if 'Per base sequence content' in eachline:
                    gc_flag = 1
                    continue
                if 'Per base N content' in eachline:
                    n_flag = 1
                    continue
                if gc_flag == 1:
                    add_gc_and_n_inf(each_sample_gc_info_dict, eachline, n)
                if n_flag == 1:
                    add_gc_and_n_inf(each_sample_gc_info_dict,
                                     eachline, n, True)
    # output each sample gc
    each_sample_gc_info_df = pd.DataFrame(each_sample_gc_info_dict)
    each_sample_gc_info_table = os.path.join(
        gc_dir, '{0}.gc.txt'.format(each_sample))
    each_sample_gc_info_df.to_csv(
        each_sample_gc_info_table, sep='\t', float_format='%.3f', columns=GC_HEADER, index=False)
    # output each sample reads quality
    with open(each_reads_quality_file, 'w') as each_reads_quality_file_inf:
        quality_list = sorted(reads_quality_dict[each_sample].keys())
        each_reads_quality_file_inf.write('Quality\tCount\tProportion\n')
        quality_count_list = []
        for each_quality in quality_list:
            quality_count_list.append(
                reads_quality_dict[each_sample][each_quality])
        for n, each_quality in enumerate(quality_list):
            each_quality_portion = float(
                quality_count_list[n]) / sum(quality_count_list)
            each_quality_count = quality_count_list[n]
            each_reads_quality_file_inf.write(
                '{each_quality}\t{each_quality_count}\t{each_quality_portion}\n'.format(**locals()))
            # merged_quality_file_inf.write(
            #     '{each_sample}\t{each_quality}\t{each_quality_count}\t{each_quality_portion}\n'.format(**locals()))
# merged_quality_file_inf.close()

out_summary_txt_file = '{}.txt'.format(out_summary_prefix)
out_summary_xlsx_file = '{}.xlsx'.format(out_summary_prefix)

header = 'Sample_ID\tReads_number(M)\tReads_length(bp)\tData_size(G)\tQ30(%)\tGC(%)'
header_list = header.split('\t')

out_summary_txt_inf = open(out_summary_txt_file, 'w')
out_summary_txt_inf.write('{}\n'.format(header))

workbook = xlsxwriter.Workbook(out_summary_xlsx_file)
worksheet = workbook.add_worksheet()
row_num = 0
output_to_xlsx(worksheet, header_list, row_num)
row_num += 1

for each_sample in sample_list:
    reads_num = round(sample_info_dict[each_sample][1] / float(10**6), 2)
    read_length = sample_info_dict[each_sample][0]
    data_size = round(sample_info_dict[each_sample][1] *
                      sample_info_dict[each_sample][0] / float(1000**3), 2)
    q30 = round(sample_info_dict[each_sample][2] *
                100 / float(sample_info_dict[each_sample][1]), 2)
    gc_content = sum(sample_info_dict[each_sample][3]) / \
        float(len(sample_info_dict[each_sample][3]))
    out_line = '%s\t%s\t%s\t%s\t%s\t%s' % (
        each_sample, reads_num, read_length, data_size, q30, gc_content)
    out_summary_txt_inf.write('{}\n'.format(out_line))
    out_list = [each_sample, reads_num,
                read_length, data_size, q30, gc_content]
    output_to_xlsx(worksheet, out_list, row_num)
    row_num += 1
out_summary_txt_inf.close()
workbook.close()
