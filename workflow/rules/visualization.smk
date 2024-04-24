rule grm_heatmap:
  input: grm = f"{config['build_directory']}/{{name}}/{{distance_method}}.grm.csv",
  output: file = f"{config['results_directory']}/{{name}}/{{distance_method}}_heatmap.png"
  script: 
    "../scripts/grm_csv_heatmap.R"


rule grm_histogram:
  input: grm = f"{config['build_directory']}/{{name}}/{{distance_method}}.grm.csv",
  output: file = f"{config['results_directory']}/{{name}}/{{distance_method}}_histogram.png"
  script: 
    "../scripts/grm_csv_histogram.R"


rule grm_to_csv:
  input: grm = f"{config['build_directory']}/{{name}}.grm",
  output: file = f"{config['build_directory']}/{{name}}.grm.csv"
  script: 
    "../scripts/grm_to_csv.R"


# snakemake --forceall --dag | dot -Tpdf > /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/GREML_dag.pdf
# snakemake --report /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/report.html
