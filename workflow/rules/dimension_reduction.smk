# rule generate_pca_variances:
#     input: 
#         data_input=config['dataset']['annotation_file'],
#     output:
#         csv_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/pca_variances.csv",
#         png_output=f"{config['results_directory']}/{config['dataset']['workname']}/pca_variances.png",
#     conda: '../../envs/dimension_reduction.yaml'
#     script: '../scripts/pca_variances.py'


# rule generate_kernelpca_variances:
#     input: 
#         data_input=config['dataset']['annotation_file'],
#     output:
#         csv_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/kernelpca_variances.csv",
#         png_output=f"{config['results_directory']}/{config['dataset']['workname']}/kernelpca_variances.png",
#     conda: '../../envs/dimension_reduction.yaml'
#     script: '../scripts/kernelpca_variances.py'


rule run_PCA:
    input: 
        input_file=config['dataset']['annotation_file'],
        pca_variances=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/pca_variances.csv",
    output:
        pca_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/{{profile,[^/]+}}.h5ad",
        nms_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/{{profile,[^/]+}}.nms",
    params:
        variance=lambda wildcards: config['profiles'][wildcards.profile]['variance'] if 'variance' in config['profiles'][wildcards.profile] else None,
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'] if 'dimensions' in config['profiles'][wildcards.profile] else None,
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
        scale="True",
    priority: 2
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/PCA.py'


rule run_kernelPCA:
    input: 
        input_file=config['dataset']['annotation_file'],
        kernelpca_variances=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/kernelpca_variances.csv",
    output:
        pca_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/{{profile,[^/]+}}.h5ad",
        nms_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/{{profile,[^/]+}}.nms",
    params:
        variance=lambda wildcards: config['profiles'][wildcards.profile]['variance'] if 'variance' in config['profiles'][wildcards.profile] else None,
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'] if 'dimensions' in config['profiles'][wildcards.profile] else None,
        kernel=lambda wildcards: config['profiles'][wildcards.profile]['kernel'],
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
        scale="True",
    priority: 2
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/kernelPCA.py'


rule grm_dim_reduction:
    input: 
        genome_input=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.method}/{wildcards.profile}.h5ad",
        nms_input=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.method}/{wildcards.profile}.nms",
    output:
        grm_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{method,pca|kernelpca}}/{{profile,[^/]+}}.grm.gz",
        grmId_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{method,pca|kernelpca}}/{{profile,[^/]+}}.grm.id",
    params:
        distance_method=lambda wildcards: config['profiles'][wildcards.profile]['distance_method'],
    priority: 1
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/generate_grm.py'
