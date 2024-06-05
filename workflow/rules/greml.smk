rule run_greml:
    input:
        grm_file=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['profiles'][wildcards.profile]['method']}/{wildcards.profile}.grm.gz",
        phenotype_file=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/phenotypes/{config['profiles'][wildcards.profile]['phenotype_file']}/{wildcards.phenotype}.phen",
    output:
        output_file=f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities/{{profile}}.{{phenotype}}.hsq",
    params:
        grm_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['profiles'][wildcards.profile]['method']}/{wildcards.profile}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities/{wildcards.profile}.{wildcards.phenotype}",
        output_folder=f"{config['build_directory']}/{config['dataset']['workname']}/heritabilities",
        maxit=250,
    wildcard_constraints:
        profile="|".join(config['profiles'].keys()),
        phenotype="|".join(config['phenotypes'].keys()),
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    shell:
        """
        mkdir -p {params.output_folder}
        if ! gcta64 --grm-gz {params.grm_prefix} --pheno {input.phenotype_file} --reml --reml-no-constrain --reml-maxit {params.maxit} --out {params.output_prefix} --thread-num {threads}; then
            echo "Error occurred, creating placeholder file."
            touch {output.output_file}
        fi
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
        file=f"{config['results_directory']}/{config['dataset']['workname']}/{config['analysis']}_heritability_summary.tsv",
    params:
        methods=[config['profiles'][profile]['method'] for profile in config['profiles']],
    conda: "../../envs/r_basic.yaml"
    script: "../scripts/heritability_summary.R"
