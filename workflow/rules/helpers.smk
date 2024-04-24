rule grm_gz_to_grm:
  input:
    grm = "{name}.grm.gz",
  output:
    grm = "{name}.grm",
  shell:
    """
    gzip -dc {input.grm} > {output.grm}
    """


rule grm_to_grm_gz:
    input:
        grm = "{file}.grm",
    output:
        grm = "{file}.grm.gz",
    shell:
        """
        gzip -c {input.grm}
        """
