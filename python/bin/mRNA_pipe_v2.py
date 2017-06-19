#! /usr/bin/python

import luigi
from os import path
from RNAseq_lib import run_cmd
from RNAseq_lib import sepcies_annotation_path
from RNAseq_lib import check_rseqc_condition
from python_tools import circ_mkdir_unix
import fastqc_pipe_v2
import quant_pipe_v2 as quant_pipe
import enrich_pipe_v2 as enrich_pipe
import star_mapping_pipe_v2
import rseqc_pipe
import snp_pipe
import rmats_pipe
import sys

script_path = path.dirname(path.abspath(__file__))


class cp_analysis_result(luigi.Task):

    from_dir = luigi.Parameter()
    to_dir = luigi.Parameter()

    def run(self):

        cp_result_cmd = ['rsync',
                         '-av',
                         '--copy-links',
                         '--exclude-from={}'.format(
                             path.join(self.from_dir, '.ignore')),
                         path.abspath(self.from_dir),
                         path.abspath(self.to_dir)]

        cp_result_inf = run_cmd(cp_result_cmd)
        with self.output().open('w') as cp_result_inf_log:
            cp_result_inf_log.write(cp_result_inf)

    def output(self):
        from_dir_name = path.basename(path.abspath(self.from_dir))
        return luigi.LocalTarget('{0}/cp_{1}_analysis_results.log'.format(log_dir, from_dir_name))


class fastqc(luigi.Task):

    dir_name = 'fastqc'

    def requires(self):
        out_dir = path.join(proj_dir, self.dir_name)
        return fastqc_pipe_v2.fastqc_collection(OutDir=out_dir, SampleInf=sample_inf, CleanDir=clean_dir)

    def run(self):
        with self.output().open('w') as fastqc_log_inf:
            fastqc_log_inf.write('fastqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/fastqc.log'.format(log_dir))


class quant(luigi.Task):

    dir_name = 'quantification'

    def requires(self):
        out_dir = path.join(proj_dir, self.dir_name)
        return quant_pipe.quant_collection(OutDir=out_dir, SampleInf=sample_inf, CleanDir=clean_dir, Transcript=transcript, Gene2Tr=gene_tr)

    def run(self):
        with self.output().open('w') as quant_log_inf:
            quant_log_inf.write('quantification finished!')

    def output(self):
        return luigi.LocalTarget('{}/quant.log'.format(log_dir))


class mapping(luigi.Task):

    dir_name = 'mapping'

    def requires(self):
        out_dir = path.join(proj_dir, self.dir_name)
        return star_mapping_pipe_v2.star_mapping_collection(OutDir=out_dir, IndexDir=star_index, SampleInf=sample_inf, CleanDir=clean_dir)

    def run(self):
        with self.output().open('w') as mapping_log_inf:
            mapping_log_inf.write('mapping finished!')

    def output(self):
        return luigi.LocalTarget('{}/mapping.log'.format(log_dir))


class enrich(luigi.Task):

    dir_name = 'enrichment'

    def requires(self):
        return quant()

    def run(self):
        out_dir = path.join(proj_dir, self.dir_name)
        quant_dir = path.join(proj_dir, 'quantification')
        yield enrich_pipe.enrichment_collection(QuantDir=quant_dir, OutDir=out_dir, GoseqAnno=goseq_ano, TopgoAnno=topgo_ano, GeneLen=gene_len, KEGGAbbr=kegg_abbr, KEGGBackground=kegg_bg, KEGGBlast=kegg_blast, ReRun='no')
        with self.output().open('w') as enrich_log_inf:
            enrich_log_inf.write('enrichment finished!')

    def output(self):
        return luigi.LocalTarget('{}/enrich.log'.format(log_dir))


class rseqc(luigi.Task):

    dir_name = 'rseqc'

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(proj_dir, 'mapping', 'bam_dir')
        out_dir = path.join(proj_dir, self.dir_name)
        yield rseqc_pipe.rseqc_collection(OutDir=out_dir, SampleInf=sample_inf, BamDir=bam_dir, BedFile=bedfile)
        with self.output().open('w') as rseqc_log_inf:
            rseqc_log_inf.write('rseqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/rseqc.log'.format(log_dir))


class snp(luigi.Task):

    dir_name = 'rseqc'

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(proj_dir, 'mapping', 'bam_dir')
        out_dir = path.join(proj_dir, self.dir_name)
        yield snp_pipe.snp_collection(OutDir=out_dir, SampleInf=sample_inf, BamDir=bam_dir, Ref=genome_fa)

        with self.output().open('w') as snp_log_inf:
            snp_log_inf.write('snp finished!')

    def output(self):
        return luigi.LocalTarget('{}/snp.log'.format(log_dir))


class splicing(luigi.Task):

    dir_name = 'splicing'

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(proj_dir, 'mapping', 'bam_dir')
        out_dir = path.join(proj_dir, self.dir_name)
        yield rmats_pipe.rmats_collection(OutDir=out_dir, SampleInf=sample_inf, BamDir=bam_dir, Gtf=gtf)

        with self.output().open('w') as splicing_log_inf:
            splicing_log_inf.write('splicing finished!')

    def output(self):
        return luigi.LocalTarget('{}/splicing.log'.format(log_dir))


class release_analysis_data(luigi.Task):

    def requires(self):

        return mapping()

    def run(self):

        analysis_bam_dir = path.join(proj_dir, 'mapping', 'bam_dir')
        out_data_dir = path.join(
            proj_dir, '{}_analysis_data'.format(proj_name))
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


class pdf_report_data(luigi.Task):

    from_dir = luigi.Parameter()
    to_dir = luigi.Parameter()

    def run(self):
        from_dir_name = path.basename(self.from_dir)
        to_dir = path.join(self.to_dir, from_dir_name)
        report_files_ini = path.join(self.from_dir, '.report_files')
        if not path.exists(report_files_ini):
            cp_cmd_inf = 'nothing for report in {}'.format(from_dir_name)
        else:
            circ_mkdir_unix(to_dir)
            cp_cmd = ['rsync',
                      '-av',
                      '--files-from={}'.format(report_files_ini),
                      self.from_dir,
                      to_dir]
            cp_cmd_inf = run_cmd(cp_cmd)
        with self.output().open('w') as cp_cmd_log:
            cp_cmd_log.write(cp_cmd_inf)

    def output(self):
        from_dir_name = path.basename(self.from_dir)
        return luigi.LocalTarget('{0}/{1}_report_data_cp.log'.format(log_dir, from_dir_name))


MODULE_DICT = {
    'fastqc': fastqc,
    'quantification': quant,
    'mapping': mapping,
    'enrichment': enrich,
    'rseqc': rseqc,
    'snp': snp,
    'splicing': splicing
}


def get_analysis_modules(module_name_list):
    module_list = []
    module_dirs = []
    for each_name in module_name_list:
        each_module = MODULE_DICT[each_name]()
        each_module_dir = path.join(proj_dir, each_name)
        module_list.append(each_module)
        module_dirs.append(each_module_dir)
    return module_list, module_dirs


class run_pipe(luigi.Task):

    # read parameter
    proj_name = luigi.Parameter()
    proj_dir = luigi.Parameter()
    clean_dir = luigi.Parameter()
    sample_inf = luigi.Parameter()
    analysis_abbr = luigi.Parameter(default='basic')
    analysis_file = luigi.Parameter(default="")
    species = luigi.Parameter()
    kegg_bg = luigi.Parameter(default='')
    database = luigi.Parameter(default='ensembl')
    database_version = luigi.Parameter(default='')

    def requires(self):

        # global paramters
        global proj_name, proj_dir, log_dir, clean_dir, sample_inf, result_dir, report_dir
        proj_name = self.proj_name
        proj_dir = self.proj_dir
        clean_dir = self.clean_dir
        sample_inf = self.sample_inf
        log_dir = path.join(proj_dir, 'logs')
        result_dir = path.join(proj_dir, '{}_analysis'.format(
            proj_name), 'analysis_result')
        report_dir = path.join(proj_dir, '{}_analysis'.format(
            proj_name), 'analysis_report')
        map(circ_mkdir_unix, [log_dir, result_dir, report_dir])

        # get species annotation
        sp_anno_inf = sepcies_annotation_path()
        sp_anno_inf.sp_latin = self.species
        sp_anno_inf.sp_database = self.database
        sp_anno_inf.sp_db_version = self.database_version
        sp_anno_inf.get_anno_inf()

        # get species annotation files
        global transcript, gene_tr, star_index, bedfile, genome_fa, gtf
        transcript = sp_anno_inf.transcript
        gene_tr = sp_anno_inf.gene_tr
        star_index = sp_anno_inf.star_index
        bedfile = sp_anno_inf.bedfile
        genome_fa = sp_anno_inf.genome_fa
        geneme_fai = sp_anno_inf.geneme_fai
        gtf = sp_anno_inf.gtf

        # enrich module annotation files
        global goseq_ano, topgo_ano, gene_len, kegg_abbr, kegg_bg, kegg_blast
        goseq_ano = sp_anno_inf.goseq_ano
        topgo_ano = sp_anno_inf.topgo_ano
        gene_len = sp_anno_inf.gene_len
        kegg_abbr = sp_anno_inf.kegg_abbr
        kegg_bg = self.kegg_bg
        if kegg_abbr == 'ko':
            if not kegg_bg:
                sys.exit('when kegg abbr is ko, kegg-bg is required!')
        else:
            kegg_bg = kegg_abbr
        kegg_blast = sp_anno_inf.kegg_blast

        # run pipeline
        global analysis_folders
        if self.analysis_file:
            analysis_list = [each.strip() for each in open(self.analysis_file)]
        else:
            if self.analysis_abbr == 'basic':
                analysis_list = ['fastqc', 'mapping',
                                 'rseqc', 'quantification', 'enrichment']
            elif self.analysis_abbr == 'advanced':
                analysis_list = ['fastqc', 'mapping',
                                 'rseqc', 'quantification',
                                 'enrichment', 'splicing', 'snp']
            else:
                sys.exit('wrong analysis_abbr!')
        ## check rseqc run condition
        if 'rseqc' in analysis_list and (not check_rseqc_condition(geneme_fai)):
            analysis_list.remove('rseqc')
        analysis_modules, analysis_folders = get_analysis_modules(
            analysis_list)
        return analysis_modules

    def run(self):
        yield [cp_analysis_result(each_folder, result_dir) for each_folder in analysis_folders]
        pdf_report_data_dir = path.join(report_dir, 'report_data')
        yield [pdf_report_data(each_folder, pdf_report_data_dir) for each_folder in analysis_folders]
        with self.output().open('w') as analysis_log:
            analysis_log.write('analysis finished')

    def output(self):
        return luigi.LocalTarget('{}/logs/analysis_finished.log'.format(self.proj_dir))


if __name__ == '__main__':
    luigi.run()
