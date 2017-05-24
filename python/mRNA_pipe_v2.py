'''

'''

import luigi
from os import path
import sys
from RNAseq_lib import run_cmd
from python_tools import circ_mkdir_unix
script_path = path.dirname(path.abspath(__file__))
fastqc_dir = path.join(script_path, 'fastqc')
quant_dir = path.join(script_path, 'quantification')
enrich_dir = path.join(script_path, 'enrichment')
star_mapping_dir = path.join(script_path, 'star_mapping')
rseqc_dir = path.join(script_path, 'rseqc')
sys.path.insert(0, fastqc_dir)
sys.path.insert(0, quant_dir)
sys.path.insert(0, enrich_dir)
sys.path.insert(0, star_mapping_dir)
sys.path.insert(0, rseqc_dir)
import fastqc_pipe_v2
import quant_pipe
import enrich_pipe
import star_mapping_pipe_v2
import rseqc_pipe


class cp_analysis_result(luigi.Task):

    from_dir = luigi.Parameter()
    to_dir = luigi.Parameter()

    def run(self):

        cp_result_cmd = ['rsync',
        '-av',
        '--copy-links',
        '--exclude-from={}'.format(path.join(self.from_dir, '.ignore')),
        path.abspath(self.from_dir),
        path.abspath(self.to_dir)]

        cp_result_inf = run_cmd(cp_result_cmd)
        with self.output().open('w') as cp_result_inf_log:
            cp_result_inf_log.write(cp_result_inf)

    def output(self):
        from_dir_name = path.basename(path.abspath(self.from_dir))
        return luigi.LocalTarget('{0}/cp_{1}_results.log'.format(log_dir, from_dir_name))


class fastqc(luigi.Task):

    def requires(self):
        return fastqc_pipe_v2.fastqc_collection(OutDir = fastqc_dir, SampleInf = sample_inf, CleanDir = clean_dir)

    def run(self):
        yield cp_analysis_result(from_dir = fastqc_dir, to_dir = result_dir)
        with self.output().open('w') as fastqc_log_inf:
            fastqc_log_inf.write('fastqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/fastqc.log'.format(log_dir))


class quant(luigi.Task):

    def requires(self):
        return quant_pipe.quant_collection(OutDir = quant_dir, SampleInf = sample_inf, CleanDir = clean_dir, Transcript = transcript, Gene2Tr = gene_tr)

    def run(self):
        yield cp_analysis_result(from_dir = quant_dir, to_dir = result_dir)
        with self.output().open('w') as quant_log_inf:
            quant_log_inf.write('quantification finished!')

    def output(self):
        return luigi.LocalTarget('{}/quant.log'.format(log_dir))


class mapping(luigi.Task):

    def requires(self):
        return star_mapping_pipe_v2.star_mapping_collection(OutDir = mapping_dir, IndexDir = star_index, SampleInf = sample_inf, CleanDir = clean_dir)

    def run(self):
        yield cp_analysis_result(from_dir = mapping_dir, to_dir = result_dir)
        with self.output().open('w') as mapping_log_inf:
            mapping_log_inf.write('mapping finished!')

    def output(self):
        return luigi.LocalTarget('{}/mapping.log'.format(log_dir))


class enrich(luigi.Task):

    def requires(self):
        return quant()

    def run(self):
        yield enrich_pipe.enrichment_collection(QuantDir = quant_dir, OutDir = enrich_dir, GoseqAnno = goseq_ano, TopgoAnno = topgo_ano, GeneLen = gene_len, KEGGAbbr = kegg_abbr, KEGGBlast = kegg_blast, ReRun = 'no')
        yield cp_analysis_result(from_dir = enrich_dir, to_dir = result_dir)
        with self.output().open('w') as enrich_log_inf:
            enrich_log_inf.write('enrichment finished!')

    def output(self):
        return luigi.LocalTarget('{}/enrich.log'.format(log_dir))


class rseqc(luigi.Task):

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(mapping_dir, 'bam_dir')
        yield rseqc_pipe.rseqc_collection(OutDir = rseqc_dir, SampleInf = sample_inf, BamDir = bam_dir, BedFile = bedfile)
        yield cp_analysis_result(from_dir = rseqc_dir, to_dir = result_dir)
        with self.output().open('w') as rseqc_log_inf:
            rseqc_log_inf.write('rseqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/rseqc.log'.format(log_dir))


class release_analysis_data(luigi.Task):

    def requires(self):

        return mapping()

    def run(self):

        analysis_bam_dir = path.join(mapping_dir, 'bam_dir')
        out_data_dir = path.join(proj_dir, 'analysis_data')
        out_bam_dir = path.join(out_data_dir, 'bam')
        fq_dir = path.join(out_data_dir, 'fq')

        circ_mkdir_unix(out_data_dir)

        ln_fq_cmd = ['ln',
        '-s',
        clean_dir,
        fq_dir]

        ln_bam_cmd = ['ln',
        '-s',
        analysis_bam_dir,
        out_bam_dir]

        link_cmd_inf = run_cmd([ln_fq_cmd, ln_bam_cmd])

        with self.output().open('w') as get_analysis_data_log:
            get_analysis_data_log.write(link_cmd_inf)

    def output(self):
        return luigi.LocalTarget('{}/release_analysis_data.log'.format(log_dir))

class snp(luigi.Task):

    def requires(self):
        return mapping()


class splicing(luigi.Task):

    def requires(self):
        return mapping()


class run_pipe(luigi.Task):

    ## read parameter
    proj_name = luigi.Parameter()
    proj_dir = luigi.Parameter()
    clean_dir = luigi.Parameter()
    sample_inf = luigi.Parameter()
    analysis = luigi.Parameter()
    transcript = luigi.Parameter()
    gene_tr = luigi.Parameter()
    goseq_ano = luigi.Parameter()
    topgo_ano = luigi.Parameter()
    gene_len = luigi.Parameter()
    kegg_abbr = luigi.Parameter()
    kegg_blast = luigi.Parameter()
    star_index = luigi.Parameter()
    bedfile = luigi.Parameter()

    def requires(self):

        ## global paramters
        global proj_name, proj_dir, log_dir, clean_dir, sample_inf, result_dir
        proj_name = self.proj_name
        proj_dir = self.proj_dir
        clean_dir = self.clean_dir
        sample_inf = self.sample_inf
        log_dir = path.join(proj_dir, 'logs')
        result_dir = path.join(proj_dir, '{}_analysis'.format(proj_name))
        map(circ_mkdir_unix, [log_dir, result_dir])

        ## fastqc module
        global fastqc_dir
        fastqc_dir = path.join(proj_dir, 'fastqc')

        ## quant module
        global quant_dir, transcript, gene_tr
        quant_dir = path.join(proj_dir, 'quantification')
        transcript = self.transcript
        gene_tr = self.gene_tr

        ## enrich module
        global enrich_dir, goseq_ano, topgo_ano, gene_len, kegg_abbr, kegg_blast
        enrich_dir = path.join(proj_dir, 'enrichment')
        goseq_ano = self.goseq_ano
        topgo_ano = self.topgo_ano
        gene_len = self.gene_len
        kegg_abbr = self.kegg_abbr
        kegg_blast = self.kegg_blast

        ## star mapping module
        global mapping_dir, star_index
        mapping_dir = path.join(proj_dir, 'mapping')
        star_index = self.star_index

        ## rseqc module
        global rseqc_dir, bedfile
        rseqc_dir = path.join(proj_dir, 'rseqc')
        bedfile = self.bedfile

        ## run pipeline
        if self.analysis == 'basic':
            return [fastqc(), rseqc(), enrich(), release_analysis_data()]
        elif self.analysis == 'advanced':
            return [fastqc(), rseqc(), enrich(), snp(), splicing(), release_analysis_data()]

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
