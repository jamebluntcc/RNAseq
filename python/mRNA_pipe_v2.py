#! /usr/bin/python

import luigi
from os import path
import sys
from RNAseq_lib import run_cmd
from RNAseq_lib import sepcies_annotation_path
from python_tools import circ_mkdir_unix
from python_tools import load_fn_to_obj

script_path = path.dirname(path.abspath(__file__))
fastqc_dir = path.join(script_path, 'fastqc')
quant_dir = path.join(script_path, 'quantification')
enrich_dir = path.join(script_path, 'enrichment')
star_mapping_dir = path.join(script_path, 'star_mapping')
rseqc_dir = path.join(script_path, 'rseqc')
snp_dir = path.join(script_path, 'snp')
splicing_dir = path.join(script_path, 'splicing')

sys.path.insert(0, fastqc_dir)
sys.path.insert(0, quant_dir)
sys.path.insert(0, enrich_dir)
sys.path.insert(0, star_mapping_dir)
sys.path.insert(0, rseqc_dir)
sys.path.insert(0, snp_dir)
sys.path.insert(0, splicing_dir)

import fastqc_pipe_v2
import quant_pipe_v2 as quant_pipe
import enrich_pipe_v2 as enrich_pipe
import star_mapping_pipe_v2
import rseqc_pipe
import snp_pipe
import rmats_pipe


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

    def requires(self):
        return fastqc_pipe_v2.fastqc_collection(OutDir=fastqc_dir, SampleInf=sample_inf, CleanDir=clean_dir)

    def run(self):
        with self.output().open('w') as fastqc_log_inf:
            fastqc_log_inf.write('fastqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/fastqc.log'.format(log_dir))


class quant(luigi.Task):

    def requires(self):
        return quant_pipe.quant_collection(OutDir=quant_dir, SampleInf=sample_inf, CleanDir=clean_dir, Transcript=transcript, Gene2Tr=gene_tr)

    def run(self):
        with self.output().open('w') as quant_log_inf:
            quant_log_inf.write('quantification finished!')

    def output(self):
        return luigi.LocalTarget('{}/quant.log'.format(log_dir))


class mapping(luigi.Task):

    def requires(self):
        return star_mapping_pipe_v2.star_mapping_collection(OutDir=mapping_dir, IndexDir=star_index, SampleInf=sample_inf, CleanDir=clean_dir)

    def run(self):
        with self.output().open('w') as mapping_log_inf:
            mapping_log_inf.write('mapping finished!')

    def output(self):
        return luigi.LocalTarget('{}/mapping.log'.format(log_dir))


class enrich(luigi.Task):

    def requires(self):
        return quant()

    def run(self):
        yield enrich_pipe.enrichment_collection(QuantDir=quant_dir, OutDir=enrich_dir, GoseqAnno=goseq_ano, TopgoAnno=topgo_ano, GeneLen=gene_len, KEGGAbbr=kegg_abbr, KEGGBlast=kegg_blast, ReRun='no')
        with self.output().open('w') as enrich_log_inf:
            enrich_log_inf.write('enrichment finished!')

    def output(self):
        return luigi.LocalTarget('{}/enrich.log'.format(log_dir))


class rseqc(luigi.Task):

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(mapping_dir, 'bam_dir')
        yield rseqc_pipe.rseqc_collection(OutDir=rseqc_dir, SampleInf=sample_inf, BamDir=bam_dir, BedFile=bedfile)
        with self.output().open('w') as rseqc_log_inf:
            rseqc_log_inf.write('rseqc finished!')

    def output(self):
        return luigi.LocalTarget('{}/rseqc.log'.format(log_dir))


class snp(luigi.Task):

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(mapping_dir, 'bam_dir')
        yield snp_pipe.snp_collection(OutDir=snp_dir, SampleInf=sample_inf, BamDir=bam_dir, Ref=genome_fa)

        with self.output().open('w') as snp_log_inf:
            snp_log_inf.write('snp finished!')

    def output(self):
        return luigi.LocalTarget('{}/snp.log'.format(log_dir))


class splicing(luigi.Task):

    def requires(self):
        return mapping()

    def run(self):
        bam_dir = path.join(mapping_dir, 'bam_dir')
        yield rmats_pipe.rmats_collection(OutDir=snp_dir, SampleInf=sample_inf, BamDir=bam_dir, Gtf=gtf)

        with self.output().open('w') as splicing_log_inf:
            splicing_log_inf.write('splicing finished!')

    def output(self):
        return luigi.LocalTarget('{}/splicing.log'.format(log_dir))


class release_analysis_data(luigi.Task):

    def requires(self):

        return mapping()

    def run(self):

        analysis_bam_dir = path.join(mapping_dir, 'bam_dir')
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
            # report_files_list = load_fn_to_obj(report_files_ini)
        #     cp_cmd_list = []
        #     for each_file in report_files_list:
        #         each_file_path = path.join(self.from_dir, each_file)
        #         cp_cmd_list.append(['cp {0} {1}'.format(each_file_path, to_dir)])
        #     cp_cmd_inf = run_cmd(cp_cmd_list, True)
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


class run_pipe(luigi.Task):

    # read parameter
    proj_name = luigi.Parameter()
    proj_dir = luigi.Parameter()
    clean_dir = luigi.Parameter()
    sample_inf = luigi.Parameter()
    analysis = luigi.Parameter(default='basic')
    species = luigi.Parameter()
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

        # fastqc module
        global fastqc_dir
        fastqc_dir = path.join(proj_dir, 'fastqc')

        # quant module
        global quant_dir, transcript, gene_tr
        quant_dir = path.join(proj_dir, 'quantification')
        transcript = sp_anno_inf.transcript
        gene_tr = sp_anno_inf.gene_tr

        # enrich module
        global enrich_dir, goseq_ano, topgo_ano, gene_len, kegg_abbr, kegg_blast
        enrich_dir = path.join(proj_dir, 'enrichment')
        goseq_ano = sp_anno_inf.goseq_ano
        topgo_ano = sp_anno_inf.topgo_ano
        gene_len = sp_anno_inf.gene_len
        kegg_abbr = sp_anno_inf.kegg_abbr
        kegg_blast = sp_anno_inf.kegg_blast

        # star mapping module
        global mapping_dir, star_index
        mapping_dir = path.join(proj_dir, 'mapping')
        star_index = sp_anno_inf.star_index

        # rseqc module
        global rseqc_dir, bedfile
        rseqc_dir = path.join(proj_dir, 'rseqc')
        bedfile = sp_anno_inf.bedfile

        # snp module
        global snp_dir, genome_fa
        snp_dir = path.join(proj_dir, 'snp')
        genome_fa = sp_anno_inf.genome_fa

        # splicing module
        global gtf
        gtf = sp_anno_inf.gtf

        # run pipeline
        global analysis_folders
        analysis_folders = []
        if self.analysis == 'basic':
            analysis_folders = [fastqc_dir, quant_dir,
                                enrich_dir, mapping_dir, rseqc_dir]
            return [fastqc(), rseqc(), enrich(), release_analysis_data()]
        elif self.analysis == 'advanced':
            analysis_folders = [fastqc_dir, quant_dir,
                                enrich_dir, mapping_dir, rseqc_dir]
            return [fastqc(), rseqc(), enrich(), snp(), release_analysis_data()]

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
