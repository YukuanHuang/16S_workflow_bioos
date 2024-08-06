version 1.0

import "./tasks/tasks_qiime2.wdl" as qiime 
import "./tasks/task_picurst2.wdl" as picurst2
import "./tasks/task_report.wdl" as report

workflow amplicon16S {
    
    input {
        Array[File] input_fastqs
        Array[File] input_fastqs_r = []
        String input_format
        String input_type
        File metadata

        String ncpus
        String qiime2_docker_image
        String picrust2_docker_image
        String memory_usage

        String p_adapter
        String p_front
        String p_adapter_r 
        String p_front_r
        String p_chimera_method
        String p_trunc_len
        String? p_trunc_len_r
        String? p_trim_left
        String? p_trunc_q

        String qiime2_min_frequency 
        String qiime2_min_samples

        File reference_sequences
        File reference_taxonomy
        String? min_length
        String? max_length
        File? pretrained_classifier

        Float? mask_max_gap_frequency
        Float? mask_min_conservation
        String parttree

        String formula 
        Int? max_nsti
        Int? min_reads
        Int? min_samples 

        String report_docker_image
        File report_assets
    }
    if(input_format == "SingleEndFastqManifestPhred33V2"){
        call qiime.qiime2_tools_import_se {
            input:
                input_fastqs = input_fastqs,
                input_format = input_format,
                input_type = input_type,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage
        }
        call qiime.qiime2_cutadapt_trim_single {
            input:
                in_data = qiime2_tools_import_se.demux_qza,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage,
                p_adapter = p_adapter,
                p_front = p_front
        }
        call qiime.qiime2_dada2_denoise_single {
            input: 
                clean_reads = qiime2_cutadapt_trim_single.trimmed_sequences,
                p_chimera_method = p_chimera_method,
                p_trunc_len = p_trunc_len,
                p_trim_left = p_trim_left,
                p_trunc_q = p_trunc_q,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage
        }
    }
    if(input_format == "PairedEndFastqManifestPhred33V2"){
        call qiime.qiime2_tools_import_pe {
            input:
                input_fastqs = input_fastqs,
                input_fastqs_r = input_fastqs_r,
                input_format = input_format,
                input_type = input_type,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage
        }
        call qiime.qiime2_cutadapt_trim_paired {
            input:
                in_data = qiime2_tools_import_pe.demux_qza,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage,
                p_adapter_f = p_adapter,
                p_front_f = p_front,
                p_adapter_r = p_adapter_r,
                p_front_r = p_front_r
        }
        call qiime.qiime2_dada2_denoise_paired {
            input: 
                clean_reads = qiime2_cutadapt_trim_paired.trimmed_sequences,
                p_chimera_method = p_chimera_method,
                p_trunc_len_f = p_trunc_len,
                p_trunc_len_r = select_first([p_trunc_len_r,p_trunc_len]),
                p_trim_left = p_trim_left,
                p_trunc_q = p_trunc_q,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage
        }
    }
    #_________________________________________
    if(!defined(pretrained_classifier)){    
        call qiime.train_classifier {
            input: 
                reference_sequences = reference_sequences,
                reference_taxonomy = reference_taxonomy,
                forward_adapter = p_front,
                reverse_adapter = p_adapter,
                min_length = min_length,
                max_length = max_length,
                p_trunc_len = p_trunc_len,
                ncpus = ncpus,
                qiime2_docker_image = qiime2_docker_image,
                memory_usage = memory_usage
        }
    }
    call qiime.qiime_filter{
        input:
            seq = select_first([qiime2_dada2_denoise_single.seqs,qiime2_dada2_denoise_paired.seqs]),
            table = select_first([qiime2_dada2_denoise_single.table,qiime2_dada2_denoise_paired.table]),
            min_frequency = qiime2_min_frequency,
            min_samples = qiime2_min_samples,
            qiime2_docker_image = qiime2_docker_image
        
    }
    call qiime.tax_analysis {
        input:
            trained_classifier = select_first([pretrained_classifier, train_classifier.trained_classifier]),
            representative_seqs_qza = qiime_filter.filtered_seqs,
            representative_table_qza = qiime_filter.filtered_table,
            ncpus = ncpus,
            qiime2_docker_image = qiime2_docker_image,
            memory_usage = memory_usage
    }
    call qiime.align_to_tree_mafft_fasttree {
        input:
            sequences = qiime_filter.filtered_seqs,
            mask_max_gap_frequency = mask_max_gap_frequency,
            mask_min_conservation = mask_min_conservation,
            parttree = parttree,
            ncpus = ncpus,
            qiime2_docker_image = qiime2_docker_image,
            memory_usage = memory_usage
    }
    call qiime.qiime2_tabulate_min_max{
        input:
            table = qiime_filter.filtered_table,
            qiime2_docker_image = qiime2_docker_image
    }
    call qiime.qiime2_diversity_core{
        input:
            metadata = metadata,
            phylogeny_tree = align_to_tree_mafft_fasttree.rooted_tree,
            table = qiime_filter.filtered_table,
            ncpus = ncpus,
            qiime2_docker_image = qiime2_docker_image,
            memory_usage = memory_usage
    }
    call qiime.qiime2_ancombc{
        input:
            metadata = metadata,
            table = qiime_filter.filtered_table,
            formula = formula,
            taxonomy = tax_analysis.tax_classification,
            qiime2_docker_image = qiime2_docker_image
    }
    call qiime.qiime2_export_for_picurst2{
        input:
            representative_seqs_qza = qiime_filter.filtered_seqs,
            representative_table_qza = qiime_filter.filtered_table,
            qiime2_docker_image = qiime2_docker_image
    }
    call picurst2.pircust2_pipeline{
        input:
            study_seqs_fna = qiime2_export_for_picurst2.rep_fasta,
            study_seqs_biom = qiime2_export_for_picurst2.feature_table,
            max_nsti = max_nsti,
            min_reads = min_reads,
            min_samples = min_samples,
            picrust2_docker_image = picrust2_docker_image,
            ncpus = ncpus,
            memory_usage = memory_usage
    }
    call report.quarto_render{
        input:
            qiime2_pcoa = qiime2_diversity_core.pcoa,
            qiime2_faith_pd_vector = qiime2_diversity_core.faith_pd_vector,
            qiime2_evenness_vector = qiime2_diversity_core.evenness_vector,
            qiime2_shannon_vector = qiime2_diversity_core.shannon_vector,
            qiime2_dada_deno_stat = select_first([qiime2_dada2_denoise_single.denoised,qiime2_dada2_denoise_paired.denoised]),
            table = qiime_filter.filtered_table,
            seq = qiime_filter.filtered_seqs,
            taxonomy = tax_analysis.tax_classification,
            tree = align_to_tree_mafft_fasttree.rooted_tree,
            picrust2_metacyc = pircust2_pipeline.metacyc
    }
}