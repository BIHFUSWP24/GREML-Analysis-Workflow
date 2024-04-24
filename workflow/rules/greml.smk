rule run_greml:
    input:
        grm_file=f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.grm.gz",
        phenotype_file=f"{config['build_directory']}/phenotypes/{{phenotype}}.phen",
    output:
        output_file=f"{config['build_directory']}/{{data_set}}_f/{{phenotype}}.hsq",
    params:
        grm_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}_f/{wildcards.data_set}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}_f/{wildcards.phenotype}",
        results_directory=config['results_directory'],
    wildcard_constraints:
        data_set="[^/]+",
        phenotype="[^/]+",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    shell:
        """
        mkdir -p {params.results_directory}/{wildcards.data_set}_f
        gcta64 --grm-gz {params.grm_prefix} --pheno {input.phenotype_file} --reml --out {params.output_prefix} --thread-num {threads}
        """


rule summerize_heritabilities:
    input:
        hsq=expand("{build_directory}/{{data_set}}_f/{phenotype}_{mode}.hsq", build_directory=config['build_directory'], phenotype=config['phenotypes'], mode=["original", "normalized"]),
    output: 
        file=f"{config['results_directory']}/{{data_set}}_greml_summary.tsv",
    params:
        build_directory=config['build_directory'],
        phenotypes=config['phenotypes'],
    conda: "../../envs/r_basic.yaml"
    script: "../scripts/heritability_summary.R"
