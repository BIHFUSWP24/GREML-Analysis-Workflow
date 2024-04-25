rule filter_genome:
    input:
        bgen=lambda wildcards: f"{config['dataset']['directory']}/{config['dataset']['format']}.bgen",
        sample=lambda wildcards: f"{config['dataset']['directory']}/{config['dataset']['format']}.sample",
        keep_fam=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['dataset']['keep_fam']['file']}",
        extract=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{config['dataset']['extract']['file']}",
    output:
        output_bed=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.bed",
        output_bim=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.bim",
        output_fam=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.fam",
    params:
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{wildcards.chromosome}",
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
        input_bed=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.bed",
        input_bim=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.bim",
        input_fam=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome}}.fam",
    output:
        output_grm_bin=f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{{chromosome}}.grm.gz",
        output_grm_id=f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{{chromosome}}.grm.id",
    params:
        input_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{wildcards.chromosome}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{wildcards.chromosome}",
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


rule merge_grm:
    input: 
        grm=lambda wildcards: expand("{directory}/{workname}/grm/chromosomes/chr{chromosome}.grm.gz",
            directory=config['build_directory'],
            workname=config['dataset']['workname'],
            chromosome=config['profiles'][wildcards.profile]['chromosomes']),
        grmId=lambda wildcards: expand("{directory}/{workname}/grm/chromosomes/chr{chromosome}.grm.id",
            directory=config['build_directory'],
            workname=config['dataset']['workname'],
            chromosome=config['profiles'][wildcards.profile]['chromosomes']),
        selection=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/multi-grm-selections/{wildcards.profile}.txt",
    output:
        grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{profile}}.grm.gz",
        grmId=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{profile}}.grm.id",
    params:
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.profile}",
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
