# gcta64 --grm-gz test_merged --pheno /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/testData/test.phen --reml --out /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/test
rule run_greml:
  input:
    grm=f"{config['build_directory']}/{{name}}_merged.grm",
    grmId=f"{config['build_directory']}/{{name}}_merged.grm.id",
    phenotypes=f"{config['data_directory']}/{{name}}.phen",
  output: f"{config['results_directory']}/{{name}}.hsq",
  params:
    in_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.name}_merged",
    out_prefix=lambda wildcards: f"{config['results_directory']}/{wildcards.name}",
  threads: config['threads']
  shell:
    """
    gcta64 --grm-gz {params.in_prefix} --pheno {input.phenotypes} --reml --out {params.out_prefix} --thread-num {threads}
    """
