# plink2 --pmerge-list /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3_f/multi-bed-selection.txt --make-bed --out /sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3_f/ukb_imp_v3
# rule combine_datasets:
#     input:
#         input_bed=lambda wildcards: expand("{build_directory}/{workname}/bfiles/chr{chromosome}.bed", build_directory=config['build_directory'], workname=config['dataset']['workname'], chromosome=config['profiles'][wildcards.profile]['chromosomes']),
#         input_bim=lambda wildcards: expand("{build_directory}/{workname}/bfiles/chr{chromosome}.bim", build_directory=config['build_directory'], workname=config['dataset']['workname'], chromosome=config['profiles'][wildcards.profile]['chromosomes']),
#         input_fam=lambda wildcards: expand("{build_directory}/{workname}/bfiles/chr{chromosome}.fam", build_directory=config['build_directory'], workname=config['dataset']['workname'], chromosome=config['profiles'][wildcards.profile]['chromosomes']),
#         selection=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/multi-bed-selection/{wildcards.profile}.txt",
#     output:
#         output_bed=f"{config['build_directory']}/{config['dataset']['workname']}/combinedBflies/{{profile,[^/]+}}.bed",
#         output_bim=f"{config['build_directory']}/{config['dataset']['workname']}/combinedBflies/{{profile,[^/]+}}.bim",
#         output_fam=f"{config['build_directory']}/{config['dataset']['workname']}/combinedBflies/{{profile,[^/]+}}.fam",
#     conda: '../../envs/embeddings.yaml'
#     params:
#         output_prefix=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/combinedBflies/{wildcards.profile}",
#     shell:
#         """
#         plink2 --pmerge-list {input.selection} --make-bed --out {params.output_prefix}
#         """


rule run_PCA:
    input: 
        input_file=config['dataset']['annotation_file']
    output:
        output_file=f"{config['build_directory']}/{config['dataset']['workname']}/pca/{{profile,[^/]+}}.h5ad",
    params:
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'],
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/PCA.py'


rule run_kernelPCA:
    input: 
        input_file=config['dataset']['annotation_file']
    output:
        output_file=f"{config['build_directory']}/{config['dataset']['workname']}/kernelpca/{{profile,[^/]+}}.h5ad",
    params:
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'],
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/kernelPCA.py'


rule grm_dim_reduction:
    input: 
        genome_input=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/{wildcards.manipulation}/{wildcards.profile}.h5ad",
    output:
        grm=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{manipulation,pca|kernelpca}}/{{profile,[^/]+}}.grm.gz",
        grmId=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{manipulation,pca|kernelpca}}/{{profile,[^/]+}}.grm.id",
    params:
        distance_method=lambda wildcards: config['profiles'][wildcards.profile]['distance_method'],
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/generate_grm.py'
