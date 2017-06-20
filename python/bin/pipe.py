#! /usr/bin/python

'''
Usage:
    pipe

'''

from os import getcwd
from os import path
from os import system

class run_pipe:
    def __init__(self):
        self.script = ''
        self.proj_dir = ''
        self.clean_dir = ''
        self.sample_inf = ''
        self.species = ''
        self.proj_name = ''
        self.proj_ini = ''
        self.sample_number = 0

    def generate_sh(self):

        worker_number = min(self.sample_number, 10)

        script_inf = '''
#! /bin/sh
mRNA_pipe_v2.py run_pipe \\
    --proj-name {0} \\
    --proj-dir {1} \\
    --clean-dir {2} \\
    --sample-inf {3} \\
    --species {4} \\
    --workers {5} \\
    '''.format(self.proj_name, self.proj_dir, self.clean_dir, self.sample_inf, self.species, worker_number)
        if path.exists(self.proj_ini):
            script_inf = '{0}--analysis-file {1} '.format(script_inf, self.proj_ini)
        with open(self.script, 'w') as script_cont:
            script_cont.write(script_inf)

    def run_script(self):

        system('nohup sh {} > {}.log 2>&1 &'.format(self.script, self.script))

    def run(self):
        if not path.exists(self.script):
            self.generate_sh()
        else:
            self.run_script()


if __name__ == '__main__':
    my_proj = run_pipe()
    my_proj.proj_dir = getcwd()
    proj_name_tmp = path.basename(my_proj.proj_dir)
    my_proj.sample_inf = path.join(my_proj.proj_dir, 'sample.ini')
    my_proj.clean_dir = path.join(my_proj.proj_dir, 'cleandata')
    my_proj.sample_number = len(open(my_proj.sample_inf).readlines())
    my_proj.species = proj_name_tmp.split('-')[2]
    my_proj.proj_name = 'OM_mRNA_{}_{}'.format(my_proj.sample_number, my_proj.species)
    my_proj.script = path.join(my_proj.proj_dir, '{}.sh'.format(my_proj.proj_name))
    my_proj.proj_ini = path.join(my_proj.proj_dir, 'project.ini')
    my_proj.run()
