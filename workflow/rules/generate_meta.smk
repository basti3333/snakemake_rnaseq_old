import math
import platform
"""Generate all the meta data files"""

### from: https://github.com/Hoohm/dropSeqPipe/blob/master/rules/generate_meta.smk

#Which rules will be run on the host computer and not sent to nodes
localrules:
    fix_annotation_biotype,
    create_dict,
    reduce_gtf,
    create_refFlat,
    create_intervals

rule fix_annotation_biotype:
    input:
        genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    output:
        genome_dir+"/{species}_{build}_{release}/annotation_fixed_biotype.gtf"
    shell:
        """
        sed 's/gene_type/gene_biotype/g' {input} {output}
        """

rule create_dict:
    input:
        genome_dir+"/{species}_{build}_{release}/genome.fa"
    output:
        genome_dir+"/{species}_{build}_{release}/genome.dict"
    threads:1
    params:
        temp_directory=config['LOCAL']['temp_directory']
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """picard CreateSequenceDictionary -Djava.io.tmpdir={params.temp_directory} \
        REFERENCE={input}\
        OUTPUT={output}
        """

rule curate_annotation:
    input:
        biotypes=config['PARAMS']['biotypes']['gtf_biotypes'],
        annotation=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    output:
        genome_dir+"/{species}_{build}_{release}/curated_annotation.gtf"
    params:
        patterns='|'.join(config['biotypes'])
    shell:
        """cat {input.annotation} | grep -E "{params.patterns}" > {output}"""


rule reduce_gtf:
    input:
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict",
        annotation=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=temp_dir
    output:
        genome_dir+"/{species}_{build}_{release}/reduced_annotation.gtf"
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ReduceGtf -m {params.memory}\
        GTF={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        IGNORE_FUNC_TYPE='null'\
        ENHANCE_GTF='false'"""


rule create_refFlat:
    input:
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict",
        annotation=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp_directory']
    output:
        genome_dir+"/{species}_{build}_{release}/annotation.refFlat"
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ConvertToRefFlat -m {params.memory}\
        ANNOTATIONS_FILE={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule create_intervals:
    input:
        annotation_reduced=genome_dir+"/{species}_{build}_{release}/curated_reduced_annotation.gtf",
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict"
    params:
        memory=config['LOCAL']['memory'],
        reference_directory=genome_dir+"/{species}_{build}_{release}",
        temp_directory=config['LOCAL']['temp_directory'],
        prefix=genome_dir+"/{species}_{build}_{release}/annotation"
    output:
        intervals=genome_dir+"/{species}_{build}_{release}/annotation.rRNA.intervals"
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && CreateIntervalsFiles -m {params.memory}\
        REDUCED_GTF={input.annotation_reduced}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        O={params.reference_directory}\
        PREFIX={params.prefix}
        """

rule create_star_index:
    input:
        reference_file=genome_dir+"/{species}_{build}_{release}/genome.fa",
        annotation_file=genome_dir+"/{species}_{build}_{release}/annotation.gtf"  
    params:
        genomeDir=genome_dir+"/{species}_{build}_{release}/STAR_INDEX/"
    output:
        genome_dir+"/{species}_{build}_{release}//STAR_INDEX/SA"
    threads: 24
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """
        mkdir -p {params.genomeDir};\
        STAR\
        --runThreadN {threads}\
        --runMode genomeGenerate\
        --genomeDir {params.genomeDir}\
        --genomeFastaFiles {input.reference_file}\
        --sjdbGTFfile {input.annotation_file}
        """
