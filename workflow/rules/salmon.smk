rule salmon_index:
    input:
        fasta=config["reference"]["transcriptome_fasta"]
    output:
        idx_dir=directory(config["paths"]["salmon_index"])
    shell:
        """
        salmon index -t {input.fasta} -i {output.idx_dir}
        """


rule salmon_quant:
    input:
        idx_dir=rules.salmon_index.output.idx_dir,
        r1=lambda wc: SAMPLES_DF.set_index("sample").loc[wc.sample, "fastq1"],
        r2=lambda wc: SAMPLES_DF.set_index("sample").loc[wc.sample, "fastq2"]
    output:
        quant=f"{config['paths']['salmon_quant_dir']}/{{sample}}/quant.sf"
    params:
        outdir=lambda wc: f"{config['paths']['salmon_quant_dir']}/{wc.sample}",
        libtype=config["params"]["salmon_libtype"]
    shell:
        """
        mkdir -p {params.outdir}
        salmon quant \
          -i {input.idx_dir} \
          -l {params.libtype} \
          -1 {input.r1} \
          -2 {input.r2} \
          -p 4 \
          -o {params.outdir}
        """
