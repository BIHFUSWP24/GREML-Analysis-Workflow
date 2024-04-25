rule run_greml:
    input:
        grm_file=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.profile}.grm.gz",
        phenotype_file=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/phenotypes/{wildcards.phenotype}.{config['profiles'][wildcards.profile]['phenotype_file']}.phen",
    output:
        output_file=f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities/{{profile}}.{{phenotype}}.hsq",
    params:
        grm_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.profile}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities/{wildcards.profile}.{wildcards.phenotype}",
        output_folder=f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities",
    wildcard_constraints:
        data_set="[^/.]+",
        phenotype="[^/.]+",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    shell:
        """
        mkdir -p {params.output_folder}
        gcta64 --grm-gz {params.grm_prefix} --pheno {input.phenotype_file} --reml --out {params.output_prefix} --thread-num {threads}
        """


rule summerize_heritabilities:
    input:
        hsq=expand("{build_directory}/{workname}/heritabilities/{profile_phenotype}.hsq",
                   build_directory=config['build_directory'],
                   workname=config['dataset']['workname'],
                   profile_phenotype=[
                        f"{profile}.{phenotype}"
                        for profile in config['profiles']
                        for phenotype in config['profiles'][profile]['phenotypes']])
    output: 
        file=f"{config['results_directory']}/{{name}}_greml_summary.tsv",
    params:
        build_directory=config['build_directory'],
        phenotypes=config['phenotypes'],
    conda: "../../envs/r_basic.yaml"
    script: "../scripts/heritability_summary.R"
