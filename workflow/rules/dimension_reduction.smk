rule PCA:
    input: input_file = f"{config['build_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{name}}/PCA/{{dimensions}}dims.h5ad"
    script: '../scripts/embeddings/PCA.py'


rule kernelPCA:
    input: input_file = f"{config['build_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{name}}/kernelPCA/{{dimensions}}dims.h5ad"
    script: '../scripts/embeddings/kernelPCA.py'


rule tSNE:
    input: input_file = f"{config['build_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{name}}/tSNE/{{dimensions}}dims{{perplexity}}pxt.h5ad"
    script: '../scripts/embeddings/tSNE.py'


rule UMAP:
    input: input_file = f"{config['build_directory']}/{{name}}.h5ad"
    output: output_file = f"{config['build_directory']}/{{name}}/UMAP/{{dimensions}}dims{{neighbors}}nn.h5ad"
    script: '../scripts/embeddings/UMAP.py'


rule distances:
    input: input_file = f"{config['build_directory']}/{{dir}}.h5ad"
    output:
        grm_output = f"{config['build_directory']}/{{dir}}/{{distance_method}}.grm",
        gz_output = f"{config['build_directory']}/{{dir}}/{{distance_method}}.grm.gz",
    log: 'logs/{dir}/{distance_method}.log'
    script: '../scripts/distance_matrix.py'