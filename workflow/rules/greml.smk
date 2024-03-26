rule run_greml:
  input:
    grm=f"{config['directory']}/{{name}}_merged.grm.bin",
    grm_id=f"{config['directory']}/{{name}}_merged.grm.id",
    phenotypes=f"{config['directory']}/{{name}}.phen",
  output: f"{config['directory']}/{{name}}.hsq",
  params:
    prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}"
  threads: config['threads']
  shell:
    """
    gcta64 --grm {params.prefix} --pheno {input.phenotypes} --reml --out {params.prefix} --thread-num {threads}
    """
