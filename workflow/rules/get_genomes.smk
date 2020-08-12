from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

localrules:
     download_genome,
     download_annotation

rule download_genome:
    input:
        FTP.remote(config['META']['genome'][species]['genome_download_path'])
    output:
        genome_dir+"/{species}_{build}_{release}/genome.fa"
    shell:
        """gunzip -d -c {input} > {output}"""

rule download_annotation:
    input:
        FTP.remote(config['META']['genome'][species]['annotation_download_path'])
    output:
        genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    shell:
        """
        gunzip -d -c {input} > {output} &&
        """

rule index_genome:
    input:
        genome_dir+"/{species}_{build}_{release}/genome.fa"
    output:
        genome_dir+"/{species}_{build}_{release}/genome.fa.fai"
    conda:
        'envs/sm_preprocess.yaml'
    shell:
        '''
        samtools faidx {input} > {output}
        '''