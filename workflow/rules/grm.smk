# rule merge_grm:
#   input: 
#     grm = expand("{directory}/{{name}}.part_{number}_{part}.grm.bin", part=range(1, config['split_in_parts']+1), number=config['split_in_parts'], directory=config['directory']),
#     grmN = expand("{directory}/{{name}}.part_{number}_{part}.grm.N.bin", part=range(1, config['split_in_parts']+1), number=config['split_in_parts'], directory=config['directory']),
#     grmId = expand("{directory}/{{name}}.part_{number}_{part}.grm.id", part=range(1, config['split_in_parts']+1), number=config['split_in_parts'], directory=config['directory']),
#   output:
#     grm = f"{config['directory']}/{{name}}.grm.bin",
#     grmN = f"{config['directory']}/{{name}}.grm.N.bin",
#     grmId = f"{config['directory']}/{{name}}.grm.id",
#   params:
#     prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}",
#     number=config['split_in_parts'],
#   shell:
#     """
#     cat {params.prefix}.part_{params.number}_*.grm.bin > {output.grm}
#     cat {params.prefix}.part_{params.number}_*.grm.N.bin > {output.grmN}
#     cat {params.prefix}.part_{params.number}_*.grm.id > {output.grmId}
#     """


# rule grm_parts:
  # input: 
  #   f"{config['directory']}/{{name}}.bed",
  #   f"{config['directory']}/{{name}}.bim",
  #   f"{config['directory']}/{{name}}.fam",
  # output:
  #   f"{config['directory']}/{{name}}.part_{{number}}_{{part}}.grm.bin",
  #   f"{config['directory']}/{{name}}.part_{{number}}_{{part}}.grm.N.bin",
  #   f"{config['directory']}/{{name}}.part_{{number}}_{{part}}.grm.id",
  # params:
  #   prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}",
  #   autosome="--autosome" if config['use_autosomes'] else "",
  # threads: config['threads'],
  # shell:
  #   """
  #   gcta64 --bfile {params.prefix} {params.autosome} --make-grm-part {wildcards.number} {wildcards.part} --thread-num {threads} --out {params.prefix}
  #   """


rule grm_chrom:
  input: 
    f"{config['directory']}/{{name}}.bed",
    f"{config['directory']}/{{name}}.bim",
    f"{config['directory']}/{{name}}.fam",
  output:
    f"{config['directory']}/{{name}}_chr{{chromosome}}.grm.bin",
    f"{config['directory']}/{{name}}_chr{{chromosome}}.grm.N.bin",
    f"{config['directory']}/{{name}}_chr{{chromosome}}.grm.id",
  params:
    input_prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}",
    output_prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}_chr{wildcards.chromosome}",
    threads=config['threads'],
  shell:
    """
    gcta64 --bfile {params.input_prefix} --chr {wildcards.chromosome} --make-grm --thread-num {params.threads} --out {params.output_prefix}
    """


rule merge_grm:
  input: 
    grm = expand("{directory}/{{name}}_chr{chromosome}.grm.bin", chromosome=config['chromosomes'], directory=config['directory']),
    grmN = expand("{directory}/{{name}}_chr{chromosome}.grm.N.bin", chromosome=config['chromosomes'], directory=config['directory']),
    grmId = expand("{directory}/{{name}}_chr{chromosome}.grm.id", chromosome=config['chromosomes'], directory=config['directory']),
    selection = f"{config['directory']}/multi-grm-selection_{{name}}.txt",
  output:
    grm = f"{config['directory']}/{{name}}_merged.grm.bin",
    grmN = f"{config['directory']}/{{name}}_merged.grm.N.bin",
    grmId = f"{config['directory']}/{{name}}_merged.grm.id",
  params:
    prefix=lambda wildcards: f"{config['directory']}/{wildcards.name}_merged",
  shell:
    """
    gcta64 --mgrm {input.selection} --make-grm --out {params.prefix}
    """


rule write_grm_selection:
  output: f"{config['directory']}/multi-grm-selection_{{prefix}}.txt"
  params:
    chromosomes = config['chromosomes'],
    directory = config['directory'],
  shell:
    """
    prefix={wildcards.prefix}
    directory={params.directory}
    for chrom in {params.chromosomes};
    do
      echo "${{directory}}/${{prefix}}_chr${{chrom}}" >> {output}
    done
    """


rule grm_to_csv:
  input: f"{config['directory']}/{{name}}.grm.bin"
  output: file = f"{config['directory']}/{{name}}_size={{size}}.grm.csv"
  params:
    prefix = lambda wildcards: f"{config['directory']}/{wildcards.name}",
  script: 
    "../scripts/grm_io.R"

