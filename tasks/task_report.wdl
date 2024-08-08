version 1.0

task quarto_render {
    input{
        File report_assets 
        File table 
        File taxonomy
        File tree
        File seq  

        File qiime2_pcoa
        File qiime2_shannon_vector
        File qiime2_faith_pd_vector
        File qiime2_evenness_vector
        File qiime2_dada_deno_stat
        File picrust2_metacyc
        
        String ncpus
        String report_docker_image
        String memory_usage
    }

    runtime {
        docker: report_docker_image
        memory_mb: memory_usage
        cpu: ncpus       
    }

    command <<<
        unzip ~{report_assets} -d ./
        quarto render --to html
        zip -r report.zip ./_site 
    >>>

    output {
        File htmp_report_zipped = "report.zip"
    }
}
