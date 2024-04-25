rule write_bed_selection:
    output: f"{config['build_directory']}/{{data_set}}/multi-bed-selection.txt",
    params:
        chromosomes = config['chromosomes'],
        build_directory = config['build_directory'],
    shell:
        """
        data_set={wildcards.data_set}
        build_directory={params.build_directory}
        for chrom in {params.chromosomes};
        do
            echo "${{build_directory}}/${{data_set}}/chr${{chrom}}.bed ${{build_directory}}/${{data_set}}/chr${{chrom}}.bim ${{build_directory}}/${{data_set}}/chr${{chrom}}.fam" >> {output}
        done
        """


# plink2 --pmerge-list /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3_f/multi-bed-selection.txt --make-bed --out /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3_f/ukb_imp_v3
rule combine_datasets:
    input:
        input_bed = expand("{build_directory}/{{data_set}}_f/chr{chromosome}.bed", build_directory=config['build_directory'], chromosome=config['chromosomes']),
        input_bim = expand("{build_directory}/{{data_set}}_f/chr{chromosome}.bim", build_directory=config['build_directory'], chromosome=config['chromosomes']),
        input_fam = expand("{build_directory}/{{data_set}}_f/chr{chromosome}.fam", build_directory=config['build_directory'], chromosome=config['chromosomes']),
        selection = f"{config['build_directory']}/{{data_set}}_f/multi-bed-selection.txt",
    output: 
        output_bed = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.bed",
        output_bim = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.bim",
        output_fam = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}.fam",
    conda: '../../envs/embeddings.yaml'
    params:
        output_prefix = f"{config['build_directory']}/{{data_set}}_f/{{data_set}}",
        first_prefix = f"{config['build_directory']}/{{data_set}}_f/chr{config['chromosomes'][0]}",
    shell:
        """
        plink2 --pmerge-list {input.selection} --make-bed --out {params.output_prefix}
        """


rule PCA:
    input: input_file = f"{config['data_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{data_set}}_f/PCA/{{dimensions}}dims.h5ad"
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/embeddings/PCA.py'


rule kernelPCA:
    input: input_file = f"{config['build_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{name}}/kernelPCA/{{dimensions}}dims.h5ad"
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/embeddings/kernelPCA.py'


rule distances:
    input: input_file = f"{config['build_directory']}/{{dir}}.h5ad"
    output:
        grm_output = f"{config['build_directory']}/{{dir}}/{{distance_method}}.grm",
        gz_output = f"{config['build_directory']}/{{dir}}/{{distance_method}}.grm.gz",
    log: 'logs/{dir}/{distance_method}.log'
    script: '../scripts/distance_matrix.py'