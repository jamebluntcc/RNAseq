import luigi
from os import path
from os import listdir
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from RNAseq_lib import GO_ANALYSIS_R
from RNAseq_lib import KEGG_ANALYSIS_PYTHON
from RNAseq_lib import ENRICH_BARPLOT_R
# from RNAseq_lib import No_task


class prepare(luigi.Task):
    '''
    prepare directories for enrichment analysis

    '''

    OutDir = luigi.Parameter()

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
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class run_go(luigi.Task):
    '''
    run go enrichment analysis
    '''

    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['Rscript',
                       GO_ANALYSIS_R,
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

    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['python',
                       KEGG_ANALYSIS_PYTHON,
                       '--blast_out',
                       KEGGBlast,
                       '--species',
                       KEGGAbbr,
                       '--background',
                       KEGGBackground,
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
    OutDir = luigi.Parameter()

    def requires(self):
        return run_go(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['Rscript',
                       ENRICH_BARPLOT_R,
                       '--anno',
                       GoseqAnno,
                       '--table',
                       '{0}/go/{1}'.format(OutDir, self.compare),
                       '--diff',
                       '{0}/differential_analysis/{1}'.format(
                           QuantDir, self.compare),
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
    OutDir = luigi.Parameter()

    def requires(self):
        return [run_kegg(OutDir=OutDir) for each_compare in compare_list]

    def run(self):
        tmp = run_cmd(['Rscript',
                       ENRICH_BARPLOT_R,
                       '--anno',
                       KEGGBlast,
                       '--table',
                       '{0}/kegg/{1}'.format(OutDir, self.compare),
                       '--diff',
                       '{0}/differential_analysis/{1}'.format(
                           QuantDir, self.compare),
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
    KEGGBackground = luigi.Parameter(default="")
    KEGGBlast = luigi.Parameter()
    ReRun = luigi.Parameter(default="")
    # UpstreamTask = luigi.Parameter(default=No_task)

    def requires(self):
        global QuantDir, OutDir, GoseqAnno, TopgoAnno, GeneLen, KEGGAbbr, KEGGBlast, compare_list, UpstreamTask, KEGGBackground
        QuantDir = self.QuantDir
        OutDir = self.OutDir
        GoseqAnno = self.GoseqAnno
        TopgoAnno = self.TopgoAnno
        GeneLen = self.GeneLen
        KEGGAbbr = self.KEGGAbbr
        KEGGBlast = self.KEGGBlast
        if not self.KEGGBackground:
            KEGGBackground = KEGGAbbr
        else:
            KEGGBackground = self.KEGGBackground
        # UpstreamTask = self.UpstreamTask
        diff_dir = path.join(QuantDir, 'differential_analysis')
        compare_list = listdir(diff_dir)
        return [(run_kegg_barplot(compare=each_compare, OutDir=OutDir), run_go_barplot(compare=each_compare, OutDir=OutDir)) for each_compare in compare_list]

    def run(self):
        ignore_files = ['.ignore', 'logs', 'kegg/blast_out',
                        'kegg/kegg_pathway_logs', 'txt.files']
        with self.output().open('w') as ignore_inf:
            for each_file in ignore_files:
                ignore_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
