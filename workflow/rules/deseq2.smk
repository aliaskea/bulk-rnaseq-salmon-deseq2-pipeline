rule make_tx2gene:
    input:
        gtf=config["reference"]["gtf"]
    output:
        tx2gene=config["paths"]["tx2gene"]
    shell:
        """
        mkdir -p $(dirname {output.tx2gene})
        Rscript scripts/make_tx2gene.R {input.gtf} {output.tx2gene}
        """


rule deseq2:
    input:
        quant=expand(
            f"{config['paths']['salmon_quant_dir']}" + "/{sample}/quant.sf",
            sample=SAMPLES,
        ),
        samples=config["paths"]["samples"],
        contrasts=config["paths"]["contrasts"],
        tx2gene=rules.make_tx2gene.output.tx2gene
    output:
        results=f"{config['paths']['deseq2_dir']}/deseq2_results.tsv"
    params:
        deseq2_dir=config["paths"]["deseq2_dir"],
        salmon_dir=config["paths"]["salmon_quant_dir"]
    shell:
        """
        mkdir -p {params.deseq2_dir}
        Rscript scripts/run_deseq2.R \
          --samples {input.samples} \
          --contrasts {input.contrasts} \
          --tx2gene {input.tx2gene} \
          --salmon_dir {params.salmon_dir} \
          --out {output.results}
        """
