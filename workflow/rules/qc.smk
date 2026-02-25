rule fastqc:
    input:
        fastq=lambda wc: SAMPLES_DF.set_index("sample").loc[wc.sample, "fastq1"]
    output:
        html=f"{config['paths']['qc_dir']}/fastqc/{{sample}}_R1_fastqc.html",
        zip=f"{config['paths']['qc_dir']}/fastqc/{{sample}}_R1_fastqc.zip"
    params:
        qc_dir=config["paths"]["qc_dir"]
    shell:
        """
        mkdir -p {params.qc_dir}/fastqc
        fastqc {input.fastq} -o {params.qc_dir}/fastqc
        base=$(basename {input.fastq})
        base=${{base%.gz}}
        base=${{base%.fastq}}
        base=${{base%.fq}}
        mv {params.qc_dir}/fastqc/${{base}}_fastqc.html {output.html}
        mv {params.qc_dir}/fastqc/${{base}}_fastqc.zip {output.zip}
        """


rule multiqc:
    input:
        expand(
            f"{config['paths']['qc_dir']}/fastqc/{{sample}}_R1_fastqc.zip",
            sample=SAMPLES,
        )
    output:
        f"{config['paths']['qc_dir']}/multiqc_report.html"
    params:
        qc_dir=config["paths"]["qc_dir"]
    shell:
        """
        mkdir -p {params.qc_dir}
        multiqc {params.qc_dir}/fastqc -o {params.qc_dir}
        """
