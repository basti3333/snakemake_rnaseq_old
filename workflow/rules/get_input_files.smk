

rule download_ENA_fasp:
    input:
    output:
        fastq_dir+'/{sample}.fastq.gz'
    params:
        ascp="-q -QT -l 300m -P33001 -i ~/.aspera/asperaweb_id_dsa.openssh",
        fasp=lambda wildcards: get_ena_fasp(wildcards.sample),
        outDir=fastq_dir
    conda:
        'envs/sm_preprocess.yaml'
    resources: ascp_limit=1
    shell:
        '''
        ascp {params.ascp} {params.fasp} {params.outDir}
        ''' 