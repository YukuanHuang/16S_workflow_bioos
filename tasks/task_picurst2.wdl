version 1.0

task pircust2_pipeline {
    input{
        File study_seqs_fna
        File study_seqs_biom
        Int? max_nsti
        Int? min_reads
        Int? min_samples

        String ncpus
        String picrust2_docker_image
        String memory_usage
    }

    runtime {
        docker: picrust2_docker_image
        memory_mb: memory_usage
        cpu: ncpus       
    }

    command <<<
        picrust2_pipeline.py -s ~{study_seqs_fna} -i ~{study_seqs_biom} -o picrust2_out_pipeline -p ~{ncpus} \
            ~{"--max_nsti " + max_nsti } \
            ~{"--min_reads " + min_reads } \
            ~{"--min_samples " + min_samples } \
            --verbose
    >>>

    output {
        File metacyc = "picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz"
        File KO_out = "picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
        File EC_out = "picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
    }
}
