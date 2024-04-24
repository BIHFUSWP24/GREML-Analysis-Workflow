rule select_phenotype:
    input:
        input_file=lambda wildcards: config['phenotype_files'][wildcards.mode],
    output:
        output_file = f"{config['build_directory']}/phenotypes/{{phenotype}}_{{mode}}.phen",
    wildcard_constraints:
        mode = "normalized|original",
    conda:
        "../../envs/phenotype.yaml",
    script:
        "../scripts/select_phenotype.py"

