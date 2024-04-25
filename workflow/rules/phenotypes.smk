rule select_phenotype:
    input:
        input_file=lambda wildcards: config['phenotypes'][wildcards.phenotype][wildcards.phenotype_file],
    output:
        output_file = f"{config['build_directory']}/{config['dataset']['workname']}/phenotypes/{{phenotype}}.{{phenotype_file}}.phen",
    conda: "../../envs/phenotype.yaml",
    script: "../scripts/select_phenotype.py"
