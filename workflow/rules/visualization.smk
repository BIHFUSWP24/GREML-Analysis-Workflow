rule grm_heatmap:
    input: grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    output: file=f"{config['results_directory']}/{config['dataset']['workname']}/heatmap/{{profile}}.png",
    script: 
        "../scripts/grm_csv_heatmap.R"


rule grm_histogram:
    input: grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    output: file=f"{config['results_directory']}//{config['dataset']['workname']}/histogram/{{profile}}.png",
    script: 
        "../scripts/grm_csv_histogram.R"


rule grm_to_csv:
    input: grm=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.profile}.grm.gz",
    output: file=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    script: 
        "../scripts/grm_to_csv.R"


# snakemake --forceall --dag | dot -Tpdf > /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/GREML_dag.pdf
# snakemake --report /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/report.html
