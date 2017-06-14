'''
Usage:
    cp_data_by_name.py <dir_inf_file> <sample_list> <out_dir>

'''
from docopt import docopt
from os import path
from os import listdir
import time

arguments = docopt(__doc__, version = '1.0')

data_type = 'clean'
dir_inf_file_list = [each.strip() for each in open(arguments['<dir_inf_file>'])]
sample_list = [each.strip() for each in open(arguments['<sample_list>'])]
out_dir = arguments['<out_dir>']

ctime = time.strftime("%Y%m%d_%H%M%S", time.localtime())

def cp_sample_cmd(from_dir, to_dir, data_type = 'clean'):
    if data_type = 'clean':
        from_dir = path.join(from_dir, 'Cleandata')
    cp_cmd_list = ['echo "start cp {}"\ndate'.format(from_dir)]
    for each_sample in listdir(from_dir):
        if each_sample in sample_list:
            cp_cmd_list.append('cp {0}/{1} {2}'.format(from_dir, each_sample, to_dir))
    else:
        cp_cmd_list.append('echo no data in sample list')
    cp_cmd_list.append('echo "finished cp {}"\ndate'.format(from_dir))
    return cp_cmd_list

def cmd2file(cmd_list, script):
    with open(script) as script_inf:
        for eachline in cmd_list:
            script_inf.write('{}\n'.format(eachline))

if __name__ == '__main__':
    cp_cmd = []
    for each_dir in dir_inf_file_list:
        if 'Cleandata' not in listdir(each_dir):
            for dep2_each_dir in listdir(each_dir):
                dep2_each_dir_path = path.join(each_dir, dep2_each_dir)
                if 'Cleandata' in listdir(dep2_each_dir_path):
                    cp_cmd.extend(cp_sample_cmd(dep2_each_dir_path, out_dir))
        else:
            cp_cmd.extend(cp_sample_cmd(each_dir, out_dir))

    cp_script = 'cp_{}.sh'.format(ctime)
    cmd2file(cp_cmd, cp_script)
