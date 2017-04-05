import subprocess
import luigi
from os import path
from os import listdir
import ConfigParser

ConfigFile = '/home/public/scripts/RNAseq/python/configure.ini'
conf = ConfigParser.ConfigParser()
conf.read(ConfigFile)

r_home = conf.get('server34', 'r_home')
go_analysis_r = conf.get('server34', 'go_analysis_r')
kegg_analysis_python = conf.get('server34', 'kegg_analysis_python')
enrich_barplot_r = conf.get('server34', 'enrich_barplot_r')

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output

class prepare(luigi.Task):

    '''
    prepare directories for enrichment analysis

    '''

    def run(self):

        log_dir = path.join(OutDir, 'logs')
        go_dir = path.join(OutDir, 'go')
        kegg_dir = path.join(OutDir, 'kegg')

        tmp = run_cmd(['mkdir',
        '-p',
        log_dir,
        go_dir,
        kegg_dir])

        with self.output().open('w') as prepare_logs:
            prepare_logs.write('prepare finished')

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))

class run_go(luigi.Task):
    '''
    run go enrichment analysis
    '''

    def requires(self):
        return prepare()

    def run(self):
        tmp = run_cmd(['{}/Rscript'.format(r_home),
        go_analysis_r,
        '--quant_dir',
        QuantDir,
        '--go_anno',
        GoseqAnno,
        '--gene_length',
        GeneLen,
        '--topgo_anno',
        TopgoAnno,
        '--out_dir',
        OutDir])

        with self.output().open('w') as go_logs_inf:
            go_logs_inf.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/go.log'.format(OutDir))

class run_kegg(luigi.Task):
    '''
    run kegg enrichment analysis
    '''

    def requires(self):
        return run_go()

    def run(self):
        tmp = run_cmd(['python',
        kegg_analysis_python,
        '--blast_out',
        KEGGBlast,
        '--species',
        KEGGAbbr,
        '--diff_dir',
        '{}/differential_analysis/'.format(QuantDir),
        '--out_dir',
        '{}/kegg'.format(OutDir)])

        with self.output().open('w') as kegg_log_inf:
            kegg_log_inf.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/kegg.log'.format(OutDir))

class run_go_barplot(luigi.Task):
    '''
    go enrichment barplot
    '''

    compare = luigi.Parameter()

    def requires(self):
        return run_kegg()

    def run(self):
        tmp = run_cmd(['Rscript',
        enrich_barplot_r,
        '--anno',
        GoseqAnno,
        '--table',
        '{0}/go/{1}'.format(OutDir, self.compare),
        '--diff',
        '{0}/differential_analysis/{1}'.format(QuantDir, self.compare),
        '--type',
        'go',
        '--out',
        '{0}/go/{1}'.format(OutDir, self.compare)])

        with self.output().open('w') as go_plot_logs:
            go_plot_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/go_barplot.log'.format(OutDir))

class run_kegg_barplot(luigi.Task):
    '''
    kegg enrichment barplot
    '''

    compare = luigi.Parameter()

    def requires(self):
        return [run_go_barplot(compare = each_compare) for each_compare in compare_list]

    def run(self):
        tmp = run_cmd(['Rscript',
        enrich_barplot_r,
        '--anno',
        KEGGBlast,
        '--table',
        '{0}/kegg/{1}'.format(OutDir, self.compare),
        '--diff',
        '{0}/differential_analysis/{1}'.format(QuantDir, self.compare),
        '--type',
        'kegg',
        '--out',
        '{0}/kegg/{1}'.format(OutDir, self.compare)])

        with self.output().open('w') as kegg_plog_logs:
            kegg_plog_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/kegg_barplot.log'.format(OutDir))

class enrichment_collection(luigi.Task):

    QuantDir = luigi.Parameter()
    OutDir = luigi.Parameter()
    GoseqAnno = luigi.Parameter()
    TopgoAnno = luigi.Parameter()
    GeneLen = luigi.Parameter()
    KEGGAbbr = luigi.Parameter()
    KEGGBlast = luigi.Parameter()

    def requires(self):
        global QuantDir, OutDir, GoseqAnno, TopgoAnno, GeneLen, KEGGAbbr, KEGGBlast, compare_list
        QuantDir = self.QuantDir
        OutDir = self.OutDir
        GoseqAnno = self.GoseqAnno
        TopgoAnno = self.TopgoAnno
        GeneLen = self.GeneLen
        KEGGAbbr = self.KEGGAbbr
        KEGGBlast = self.KEGGBlast
        diff_dir = path.join(QuantDir, 'differential_analysis')
        compare_list = listdir(diff_dir)
        return [run_kegg_barplot(compare = each_compare) for each_compare in compare_list]

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
