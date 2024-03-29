rule grm_visualization:
  input: grm = f"{config['build_directory']}/{{name}}_merged.grm",
  output: file = f"{config['results_directory']}/{{name}}_{{mode}}.png"
  script: 
    "../scripts/grm_visualization.R"


rule grm_to_tsv:
  input: grm = f"{config['build_directory']}/{{name}}_merged.grm",
  output: file = f"{config['build_directory']}/{{name}}.tsv"
  script: 
    "../scripts/grm_to_csv.R"