# rule grm_gz_to_grm:
#   input:
#     grm = "{name}.grm.gz",
#   output:
#     grm = "{name}.grm",
#   shell:
#     """
#     gzip -dc {input.grm} > {output.grm}
#     """


rule grm_to_grm_gz:
    input:
        grm = "{file}.grm",
    output:
        grm = "{file}.grm.gz",
    shell:
        """
        gzip -c {input.grm}
        """


rule write_grm_selection:
    output: 
        grm_selection=f"{config['build_directory']}/{config['dataset']['workname']}/grm/multi-grm-selections/{{profile}}.txt"
    params:
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
        chromosome_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/multi-grm-selections/{wildcards.profile}/chr",
    shell:
        """
        for chrom in {params.chromosomes};
        do
            echo "${{params.chromosome_prefix}}${{chrom}}" >> {output.grm_selection}
        done
        """