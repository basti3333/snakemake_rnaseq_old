



rule bedgraph:
    input:
        bam='results/star/{sample}.Aligned.sortedByCoord.out.bam',
        idx=genome_dir+"/{species}_{build}_{release}/genome.fa.fai"
    output:
        'results/star/{species}_{build}_{release}/{sample}.bedgraph.gz'
    conda:
        'envs/sm_preprocess.yaml'
    shell:
        '''
        bedtools \
        genomecov \
        -split \
        -bg \
        -ibam {input.bam} \
        | gzip > {output}
        '''