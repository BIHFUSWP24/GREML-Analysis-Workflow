rule run_PCA:
    input: 
        input_file=config['dataset']['annotation_file'],
    output:
        pca_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/{{profile,[^/]+}}.h5ad",
        nms_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/pca/{{profile,[^/]+}}.nms",
    params:
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'],
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
    priority: 2
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/PCA.py'


rule run_kernelPCA:
    input: 
        input_file=config['dataset']['annotation_file'],
    output:
        pca_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/{{profile,[^/]+}}.h5ad",
        nms_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/kernelpca/{{profile,[^/]+}}.nms",
    params:
        dimensions=lambda wildcards: config['profiles'][wildcards.profile]['dimensions'],
        chromosomes=lambda wildcards: config['profiles'][wildcards.profile]['chromosomes'],
    priority: 2
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/kernelPCA.py'


rule grm_dim_reduction:
    input: 
        genome_input=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.manipulation}/{wildcards.profile}.h5ad",
        nms_input=lambda wildcards: f"{config['build_directory']}/{config['dataset']['workname']}/grm/{wildcards.manipulation}/{wildcards.profile}.nms",
    output:
        grm_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{manipulation,pca|kernelpca}}/{{profile,[^/]+}}.grm.gz",
        grmId_output=f"{config['build_directory']}/{config['dataset']['workname']}/grm/{{manipulation,pca|kernelpca}}/{{profile,[^/]+}}.grm.id",
    params:
        distance_method=lambda wildcards: config['profiles'][wildcards.profile]['distance_method'],
    priority: 1
    conda: '../../envs/dimension_reduction.yaml'
    script: '../scripts/generate_grm.py'
