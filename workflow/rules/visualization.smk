rule grm_heatmap:
    input: grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{wildcards.profile}}.grm.csv",
    output: file=f"{config['results_directory']}/{config['dataset']['workname']}/heatmap/{{profile}}.png",
    conda: "../../envs/r_plot.yaml"
    script: "../scripts/grm_csv_heatmap.R"


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
    input: summary_input=lambda wildcards: f"{config['results_directory']}/{config['dataset']['workname']}/{wildcards.name}_greml_summary.tsv",
    output: summary_output=f"{config['results_directory']}/{config['dataset']['workname']}/{{name}}_greml_summary.pdf",
    conda: "../../envs/r_plot.yaml"
    script: "../scripts/summary_plot.R"


# snakemake --forceall --dag | dot -Tpdf > /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/GREML_dag.pdf
# snakemake --report /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/report.html
