from Bio import SeqIO
import sys
import os

if not len(sys.argv) == 4 :
    print "python\t"+sys.argv[0]+"\tinput_fasta\tnumber_of_subfiles\toutput_dir\n"
    sys.exit(0)

input_fasta = sys.argv[1]
number_of_subfiles = int(sys.argv[2])
output_dir = sys.argv[3]
if not os.path.exists(output_dir):
    os.system('mkdir -p %s ' % output_dir)

seq_list = [[] for i in range(number_of_subfiles)]

for seq_num, seq_record in enumerate(SeqIO.parse(input_fasta, "fasta")) :
    number = seq_num % number_of_subfiles
    seq_list[number].append(seq_record)

for sub_num in range(number_of_subfiles):
    file_name = os.path.join(output_dir, 'split_fa_{}.fa'.format(sub_num))
    SeqIO.write(seq_list[sub_num], file_name, "fasta")
