rule get_eids_srids:
    input:
        input_file=config['dataset']['annotation_file'],
    output:
        eid_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/eid.txt",
        srid_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/srid.txt",
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/get_eids_srids.py'


rule filter_genome:
    input:
        bgen=lambda wildcards: f"{config['dataset']['directory']}/{config['dataset']['format']}.bgen",
        sample=lambda wildcards: f"{config['dataset']['directory']}/{config['dataset']['format']}.sample",
        keep_fam=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/eid.txt",
        extract=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/srid.txt",
    output:
        output_bed=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome,[0-9]+}}.bed",
        output_bim=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome,[0-9]+}}.bim",
        output_fam=f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{{chromosome,[0-9]+}}.fam",
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
        output_grm_bin=f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{{chromosome,[0-9]+}}.grm.gz",
        output_grm_id=f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{{chromosome,[0-9]+}}.grm.id",
    params:
        input_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/bfiles/chr{wildcards.chromosome}",
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/chromosomes/chr{wildcards.chromosome}",
    threads: config['threads']
    conda: "../../envs/embeddings.yaml"
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
        grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/none/{{profile,[^/]+}}.grm.gz",
        grmId=f"{config['build_directory']}/{config['dataset']['workname']}/grm/none/{{profile,[^/]+}}.grm.id",
    params:
        output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/none/{wildcards.profile}",
    threads: config['threads']
    priority: 1
    conda: "../../envs/embeddings.yaml"
    shell:
        """
        gcta64 \
            --mgrm-gz {input.selection} \
            --make-grm-gz \
            --thread-num {threads} \
            --out {params.output_prefix}
        """
