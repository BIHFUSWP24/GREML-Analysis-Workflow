rule grm_heatmap:
    input:
        grm=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['profiles'][wildcards.profile]['method']}/{wildcards.profile}.grm.gz",
        grmID=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['profiles'][wildcards.profile]['method']}/{wildcards.profile}.grm.id",
        phenotype_inputs=lambda wildcards: expand("{build_directory}/{workname}/phenotypes/{phenotype_file}/{phenotype}.phen", build_directory=config['build_directory'], workname=config['dataset']['workname'], phenotype_file=config['profiles'][wildcards.profile]['phenotype_file'], phenotype=config['profiles'][wildcards.profile]['phenotypes']),
    output: file=f"{config['results_directory']}/{config['dataset']['workname']}/heatmap/{{profile}}.png",
    params:
        phenotype_names=lambda wildcards: [config['phenotypes'][phenotype][config['profiles'][wildcards.profile]['phenotype_file']]['name'] for phenotype in config['profiles'][wildcards.profile]['phenotypes']],
    conda: "../../envs/r_heatmap.yaml"
    script: "../scripts/grm_heatmap.R"


rule grm_histogram:
    input: grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    output: file=f"{config['results_directory']}//{config['dataset']['workname']}/histogram/{{profile}}.png",
    conda: "../../envs/r_plot.yaml"
    script: "../scripts/grm_csv_histogram.R"


rule grm_to_csv:
    input: grm=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.profile}.grm.gz",
    output: file=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    conda: "../../envs/r_plot.yaml"
    script: "../scripts/grm_to_csv.R"


rule summary_plot:
    input:
        heritability_input=f"{config['results_directory']}/{config['dataset']['workname']}/{config['analysis']}_heritability_summary.tsv",
    output: summary_output=f"{config['results_directory']}/{config['dataset']['workname']}/{config['analysis']}_summary.pdf",
    conda: "../../envs/r_plot.yaml"
    script: "../scripts/heritability_plot.R"


# snakemake --forceall --dag | dot -Tpdf > /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/GREML_dag.pdf
# snakemake --report /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/report.html
