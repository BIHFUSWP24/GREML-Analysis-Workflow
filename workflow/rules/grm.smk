rule filter_genome:
    input:
        bgen = f"{config['data_directory']}/{{data_set}}_{{chromosome}}_{{version}}.bgen",
        sample = f"{config['data_directory']}/{{data_set}}_{{chromosome}}_{{version}}.sample",
        keep_fam = f"{config['build_directory']}/100k_all_eid.txt",
        extract = f"{config['build_directory']}/100k_all_srid.txt",
    output:
        output_bed = f"{config['build_directory']}/{{data_set}}_{{version}}_f/{{chromosome}}.bed",
        output_bim = f"{config['build_directory']}/{{data_set}}_{{version}}_f/{{chromosome}}.bim",
        output_fam = f"{config['build_directory']}/{{data_set}}_{{version}}_f/{{chromosome}}.fam",
    params:
        output_prefix = lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}_{wildcards.version}_f/{wildcards.chromosome}",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    priority: 2
    shell:
        """
        plink2 \
            --bgen {input.bgen} ref-first \
            --sample {input.sample} \
            --keep-fam {input.keep_fam} \
            --extract {input.extract} \
            --make-bed \
            --threads {threads} \
            --out {params.output_prefix}
        """


rule chromosome_grm:
    input: 
        input_bed = f"{config['build_directory']}/{{data_set}}/{{chromosome}}.bed",
        input_bim = f"{config['build_directory']}/{{data_set}}/{{chromosome}}.bim",
        input_fam = f"{config['build_directory']}/{{data_set}}/{{chromosome}}.fam",
    output:
        output_grm_gz = f"{config['build_directory']}/{{data_set}}/{{chromosome}}.grm.gz",
        output_grm_id = f"{config['build_directory']}/{{data_set}}/{{chromosome}}.grm.id",
    params:
        input_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}/{wildcards.chromosome}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}/{wildcards.chromosome}",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    priority: 1
    shell:
        """
        gcta64 \
            --bfile {params.input_prefix} \
            --make-grm-gz \
            --thread-num {threads} \
            --out {params.output_prefix}
        """


rule write_grm_selection:
  output: f"{config['build_directory']}/{{data_set}}/multi-grm-selection.txt",
  params:
    chromosomes = config['chromosomes'],
    build_directory = config['build_directory'],
  shell:
    """
    data_set={wildcards.data_set}
    build_directory={params.build_directory}
    for chrom in {params.chromosomes};
    do
      echo "${{build_directory}}/${{data_set}}/chr${{chrom}}" >> {output}
    done
    """


rule merge_grm:
    input: 
        grm = expand("{directory}/{{data_set}}_f/chr{chromosome}.grm.gz", chromosome=config['chromosomes'], directory=config['build_directory']),
        grmId = expand("{directory}/{{data_set}}_f/chr{chromosome}.grm.id", chromosome=config['chromosomes'], directory=config['build_directory']),
        selection = f"{config['build_directory']}/{{data_set}}_f/multi-grm-selection.txt",
    output:
        grm = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.grm.gz",
        grmId = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.grm.id",
    params:
        output_prefix=lambda wildcards: f"{config['build_directory']}/{wildcards.data_set}_f/{wildcards.data_set}",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
    shell:
        """
        gcta64 \
            --mgrm-gz {input.selection} \
            --make-grm-gz \
            --thread-num {threads} \
            --out {params.output_prefix}
        """
