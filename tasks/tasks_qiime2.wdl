version 1.0

task qiime2_tools_import_se {
    input{
        Array[File] input_fastqs
        String input_format
        String input_type
        String ncpus
        String qiime2_docker_image
        String memory_usage
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus       
    }

    command <<<
        set -ex
        echo -e "sample-id\tabsolute-filepath" > manifest.tsv
        for fastq in ~{sep=' ' input_fastqs}; do
            #making new bash variable | regex: (_) -> (-)
            SAMPLENAME=$(basename $fastq | perl -lape 's/\.(fastq|fq)(\.gz)?//g' | perl -lape 's/\.(fastq|fq)(\.gz)?//g')
            #All names added to one giant file 
            echo -e "$SAMPLENAME\t$fastq" >> manifest.tsv
        done
        qiime tools import \
            --type '~{input_type}' \
            --input-path manifest.tsv \
            --output-path demux.qza \
            --input-format ~{input_format} 
    >>>

    output {
        File demux_qza = "demux.qza"
    }
}
task qiime2_tools_import_pe {
    input{
        Array[File] input_fastqs
        Array[File] input_fastqs_r
        String input_format
        String input_type
        String ncpus
        String qiime2_docker_image
        String memory_usage
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus       
    }

    command <<<
        set -ex
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
        array1=(~{sep=' ' input_fastqs})
        array2=(~{sep=' ' input_fastqs_r})
        len=${#array1[@]}
        for(( i=0; i<len; i++ )); do
            #making new bash variable | regex: (_) -> (-)
            SAMPLENAME=$(basename ${array1[$i]} | perl -lape 's/\.(fastq|fq)(\.gz)?//g' | perl -lape 's/\.(_1|_R1)//g')
            #All names added to one giant file 
            echo -e "$SAMPLENAME\t${array1[$i]}\t${array2[$i]}" >> manifest.tsv
        done
        qiime tools import \
            --type '~{input_type}' \
            --input-path manifest.tsv \
            --output-path demux.qza \
            --input-format ~{input_format} 
    >>>

    output {
        File demux_qza = "demux.qza"
    }
}
task qiime2_demux_summarize{
    input {
        File demux_qza

        String ncpus
        String qiime2_docker_image
        String memory_usage
    }

    command <<<
       qiime demux summarize \
            --i-data ~{demux_qza} \
            --p-n ~{ncpus} \
            --o-visualization demux_summarize.qzv
    >>>

    output {
        File demux_summarize = "demux_summarize.qzv"
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus      
    }
}
task qiime2_cutadapt_trim_single {
    input {
        File in_data
        String ncpus
        String qiime2_docker_image
        String memory_usage
        String p_adapter
        String p_front
    }

    command <<<
       qiime cutadapt trim-single \
            --i-demultiplexed-sequences ~{in_data} \
            --p-cores ~{ncpus} \
            --p-adapter ~{p_adapter} \
            --p-front ~{p_front} \
            --o-trimmed-sequences "adapter_trimmed.qza"
    >>>

    output {
        File trimmed_sequences = "adapter_trimmed.qza"
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus      
    }
}
task qiime2_cutadapt_trim_paired {
    input {
        File in_data
        String ncpus
        String qiime2_docker_image
        String memory_usage
        String p_adapter_f
        String p_front_f
        String p_adapter_r
        String p_front_r
    }

    command <<<
       qiime cutadapt trim-paired \
            --i-demultiplexed-sequences ~{in_data} \
            --p-cores ~{ncpus} \
            --p-adapter-f ~{p_adapter_f} \
            --p-front-f ~{p_front_f} \
            --p-adapter-r ~{p_adapter_r} \
            --p-front-r ~{p_front_r} \
            --o-trimmed-sequences "adapter_trimmed.qza"
    >>>

    output {
        File trimmed_sequences = "adapter_trimmed.qza"
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus      
    }
}
task qiime2_dada2_denoise_single {
    input {
        File clean_reads
        String ncpus
        String qiime2_docker_image
        String memory_usage
        String p_chimera_method
        String p_trunc_len
        String? p_trim_left
        String? p_trunc_q
    }

    command <<<
        qiime dada2 denoise-single \
        --i-demultiplexed-seqs ~{clean_reads} \
        ~{"--p-trim-left " + p_trim_left} \
        ~{"--p-trunc-len " + p_trunc_len} \
        --p-n-threads ~{ncpus} \
        ~{"--p-trunc-q " + p_trunc_q} \
        --p-chimera-method ~{p_chimera_method} \
        --o-table "table.qza" \
        --o-representative-sequences "rep_seq.qza" \
        --o-denoising-stats "deno_stats.qza"
    >>>

    output {
        File seqs = "rep_seq.qza"
        File table = "table.qza"
        File denoised = "deno_stats.qza"
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus 
    }
}
task qiime2_dada2_denoise_paired {
    input {
        File clean_reads
        String ncpus
        String qiime2_docker_image
        String memory_usage
        String p_chimera_method
        String p_trunc_len_f
        String p_trunc_len_r
        String? p_trim_left
        String? p_trunc_q
    }

    command <<<
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs ~{clean_reads} \
            ~{"--p-trim-left " + p_trim_left} \
            "--p-trunc-len-f" ~{p_trunc_len_f} \
            "--p-trunc-len-r" ~{p_trunc_len_r} \
            --p-n-threads ~{ncpus} \
            ~{"--p-trunc-q " + p_trunc_q} \
            --p-chimera-method ~{p_chimera_method} \
            --o-table "table.qza" \
            --o-representative-sequences "rep_seq.qza" \
            --o-denoising-stats "deno_stats.qza"
    >>>

    output {
        File seqs = "rep_seq.qza"
        File table = "table.qza"
        File denoised = "deno_stats.qza"
    }

    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus 
    }
}
task train_classifier {
    input {
        File reference_sequences
        File reference_taxonomy
        String? forward_adapter
        String? reverse_adapter
        String? min_length
        String? max_length
        String? p_trunc_len
        String ncpus
        String qiime2_docker_image
        String memory_usage
    }
    command <<<
        set -ex
        qiime feature-classifier extract-reads \
        --i-sequences "~{reference_sequences}" \
        --p-f-primer "~{forward_adapter}" \
        --p-r-primer "~{reverse_adapter}" \
        ~{"--p-min-length " + min_length} \
        ~{"--p-max-length " + max_length} \
        ~{"--p-trunc-len " + p_trunc_len} \
        --p-n-jobs ~{ncpus} \
        --o-reads "extracts.qza"

        qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "extracts.qza" \
        --i-reference-taxonomy "~{reference_taxonomy}" \
        --o-classifier "classifier.qza"
        >>>
    output {
        File   trained_classifier = "classifier.qza"
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus
    }
}
task tax_analysis {
    input {
        File    trained_classifier
        File    representative_seqs_qza
        File    representative_table_qza
        String ncpus
        String qiime2_docker_image
        String memory_usage 
    }
    command <<<
        set -ex
        qiime feature-classifier classify-sklearn \
        --i-classifier ~{trained_classifier} \
        --i-reads ~{representative_seqs_qza} \
        --o-classification "taxonomy.qza"
        
        qiime feature-table tabulate-seqs \
            --i-data ~{representative_seqs_qza} \
            --o-visualization "list_rep_seqs.qzv"
        
        qiime feature-table summarize \
            --i-table ~{representative_table_qza} \
            --o-visualization "feature_table_summary.qzv"

        qiime taxa barplot \
            --i-table ~{representative_table_qza} \
            --i-taxonomy "taxonomy.qza" \
            --o-visualization "taxa_bar_plots.qzv"
    >>>
    output {
        File tax_classification = "taxonomy.qza"
        File rep_seq_list = "list_rep_seqs.qzv" 
        File tax_classification_graph = "taxa_bar_plots.qzv"
        File feature_table_summary = "feature_table_summary.qzv"
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus      
    }
}
task align_to_tree_mafft_fasttree {
    input {
        File sequences               
        Float? mask_max_gap_frequency = 1.0  # Default: Retains all columns regardless of gap frequency
        Float? mask_min_conservation = 0.4   # Default: Requires at least one character present in 40% of sequences
        String parttree = "False"            # Default: Disabled, enable if aligning more than 1,000,000 sequences

        String ncpus
        String qiime2_docker_image
        String memory_usage
    }
    command <<<
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences ~{sequences} \
            --p-n-threads ~{ncpus} \
            ~{"--p-mask-max-gap-frequency " + mask_max_gap_frequency} \
            ~{"--p-mask-min-conservation " + mask_min_conservation} \
            ~{"--p-parttree " + parttree} \
            --o-alignment aligned-rep-seqs.qza \
            --o-masked-alignment masked-aligned-rep-seqs.qza \
            --o-tree tree.qza \
            --o-rooted-tree rooted-tree.qza
    >>>

    output {
        File alignment = "aligned-rep-seqs.qza"     # Output: Unmasked alignment
        File masked_alignment = "masked-aligned-rep-seqs.qza" # Output: Masked alignment
        File unrooted_tree = "tree.qza"   # Output: Unrooted phylogeny
        File rooted_tree = "rooted-tree.qza"       # Output: Rooted phylogeny at midpoint
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus  
    }
    meta {
        description: "Builds a phylogenetic tree using MAFFT alignment and FastTree inference."
    }
}
task qiime2_tabulate_min_max{
    input {
        File table

        String qiime2_docker_image
    }
    command <<<
        set -ex 
        qiime feature-table tabulate-sample-frequencies \
            --i-table ~{table} \
            --o-sample-frequencies sample-frequencies.qza
        qiime tools export --input-path sample-frequencies.qza \
            --output-path ./sample_frequencies
        mv ./sample_frequencies/metadata.tsv sample_frequencies.tsv
        R --no-echo --vanilla <<Rscript
            data <- read.table("sample_frequencies.tsv", header = TRUE, sep = "\t", 
                                check.names = FALSE, stringsAsFactors = FALSE)
            min_val <- min(as.numeric(gsub(",", "", data[,2])), na.rm = TRUE)
            max_val <- max(as.numeric(gsub(",", "", data[,2])), na.rm = TRUE)
            cat(min_val, "\n")
            cat(max_val, "\n")
        Rscript
    >>>
    output{
        Array[String] values = read_lines(stdout())
        String min_val = values[0]
        String max_val = values[1]
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: "2 GB"
        cpu: 1  
    }
}
task qiime2_diversity_core {
    input {
        File metadata
        File phylogeny_tree
        String sampling_depth = 10000
        File table

        String ncpus
        String qiime2_docker_image
        String memory_usage
    }
    command <<<
        qiime diversity core-metrics-phylogenetic \
            --m-metadata-file ~{metadata} \
            --i-phylogeny ~{phylogeny_tree} \
            --i-table ~{table} \
            --p-sampling-depth ~{sampling_depth} \
            --p-n-jobs-or-threads ~{ncpus} \
            --output-dir "./diversity_core/"         
    >>>
    output {
        File rarefied_table = "./diversity_core/rarefied_table.qza"
        File faith_pd_vector = "./diversity_core/faith_pd_vector.qza"
        File observed_features_vector = "./diversity_core/observed_features_vector.qza"
        File shannon_vector = "./diversity_core/shannon_vector.qza"
        File evenness_vector = "./diversity_core/evenness_vector.qza"
        File unweighted_unifrac_distance_matrix = "./diversity_core/unweighted_unifrac_distance_matrix.qza"
        File weighted_unifrac_distance_matrix = "./diversity_core/weighted_unifrac_distance_matrix.qza"
        File jaccard_distance_matrix = "./diversity_core/jaccard_distance_matrix.qza"
        File bray_curtis_distance_matrix = "./diversity_core/bray_curtis_distance_matrix.qza"
        File unweighted_unifrac_pcoa_results = "./diversity_core/unweighted_unifrac_pcoa_results.qza"
        File weighted_unifrac_pcoa_results = "./diversity_core/weighted_unifrac_pcoa_results.qza"
        File jaccard_pcoa_results = "./diversity_core/jaccard_pcoa_results.qza"
        File bray_curtis_pcoa_results = "./diversity_core/bray_curtis_pcoa_results.qza"
        File unweighted_unifrac_emperor = "./diversity_core/unweighted_unifrac_emperor.qzv"
        File weighted_unifrac_emperor = "./diversity_core/weighted_unifrac_emperor.qzv"
        File jaccard_emperor = "./diversity_core/jaccard_emperor.qzv"
        File bray_curtis_emperor = "./diversity_core/bray_curtis_emperor.qzv"
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: memory_usage
        cpu: ncpus  
    }
}
task qiime_filter{
    input{
        File table
        File seq
        String min_frequency 
        String min_samples

        String qiime2_docker_image 
    }
    command <<<
        qiime feature-table filter-features \
            --i-table ~{table} \
            --p-min-frequency ~{min_frequency} \
            --p-min-samples ~{min_samples} \
            --o-filtered-table filtered-table.qza
    
        qiime feature-table filter-seqs \
            --i-data ~{seq} \
            --i-table filtered-table.qza \
            --o-filtered-data filtered-sequences.qza
    >>>
    output{
        File filtered_table = "filtered-table.qza"
        File filtered_seqs = "filtered-sequences.qza"
    }
    runtime{
        docker: qiime2_docker_image
        memory_mb: "4 GB"
        cpu: 1  
    }
}

task qiime2_ancombc{
    input{
        File metadata 
        File table 
        File taxonomy
        String formula 
        Array[String] levels = [2,3,4,5,6]

        String qiime2_docker_image
    }
    command <<<
        set -ex
        for tax_level in ~{sep=' ' levels}; do
            qiime taxa collapse \
                --i-table ~{table} \
                --i-taxonomy ~{taxonomy} \
                --p-level $tax_level \
                --o-collapsed-table "collapsed-table-$tax_level.qza"
            qiime composition ancombc \
                --i-table "collapsed-table-$tax_level.qza" \
                --p-formula ~{formula} \
                --m-metadata-file ~{metadata} \
                --o-differentials "dataloaf-$tax_level.qza"    
        done
    >>>
    output {
        Array[File] comp = glob("dataloaf-*.qza")
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: "2 GB"
        cpu: 1 
    }
}

task qiime2_export_for_picurst2{
    input{
        File representative_seqs_qza
        File representative_table_qza
        String qiime2_docker_image
    }
    command <<<
        set -ex
        qiime tools export \
            --input-path ~{representative_seqs_qza} \
            --output-path ./exported_data
        mv ./exported_data/* ./
        qiime tools export \
            --input-path ~{representative_table_qza} \
            --output-path ./exported_data
        mv ./exported_data/* ./
    >>>
    output {
        File feature_table = "feature-table.biom"
        File rep_fasta = "dna-sequences.fasta"
    }
    runtime {
        docker: qiime2_docker_image
        memory_mb: "2 GB"
        cpu: 1 
    }
}