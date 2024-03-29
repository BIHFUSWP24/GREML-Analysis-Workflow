# gcta64 --bfile /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/testData/test --chr 1 --make-grm-gz --out test_chr1
rule grm_chrom:
  input: 
    f"{config['data_directory']}/{{name}}.bed",
    f"{config['data_directory']}/{{name}}.bim",
    f"{config['data_directory']}/{{name}}.fam",
  output:
    f"{config['build_directory']}/{{name}}_chr{{chromosome}}.grm.gz",
    f"{config['build_directory']}/{{name}}_chr{{chromosome}}.grm.id",
  params:
    input_prefix=lambda wildcards: f"{config['data_directory']}/{wildcards.name}",
    output_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.name}_chr{wildcards.chromosome}",
    threads=config['threads'],
  shell:
    """
    gcta64 --bfile {params.input_prefix} --chr {wildcards.chromosome} --make-grm-gz --thread-num {params.threads} --out {params.output_prefix}
    """


rule write_grm_selection:
  output: f"{config['build_directory']}/multi-grm-selection_{{prefix}}.txt"
  params:
    chromosomes = config['chromosomes'],
    directory = config['build_directory'],
  shell:
    """
    prefix={wildcards.prefix}
    directory={params.directory}
    for chrom in {params.chromosomes};
    do
      echo "${{directory}}/${{prefix}}_chr${{chrom}}" >> {output}
    done
    """


# gcta64 --mgrm-gz multi-grm-selection_test.txt --make-grm-gz --out test_merged
rule merge_grm:
  input: 
    grm = expand("{directory}/{{name}}_chr{chromosome}.grm", chromosome=config['chromosomes'], directory=config['build_directory']),
    grmId = expand("{directory}/{{name}}_chr{chromosome}.grm.id", chromosome=config['chromosomes'], directory=config['build_directory']),
    selection = f"{config['build_directory']}/multi-grm-selection_{{name}}.txt",
  output:
    grm = f"{config['build_directory']}/{{name}}_merged.grm.gz",
    grmId = f"{config['build_directory']}/{{name}}_merged.grm.id",
  params:
    prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.name}_merged",
  shell:
    """
    gcta64 --mgrm-gz {input.selection} --make-grm-gz --out {params.prefix}
    """



# gzip -dc test_chr1.grm.gz > test_chr1.grm
# gzip -dc test_merged.grm.gz > test_merged.grm
rule grmgz_to_grm:
  input:
    grm = f"{config['build_directory']}/{{name}}.grm.gz",
  output:
    grm = f"{config['build_directory']}/{{name}}.grm",
  shell:
    """
    gzip -dc {input.grm} > {output.grm}
    """
