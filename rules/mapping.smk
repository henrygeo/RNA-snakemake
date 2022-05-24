import pandas as pd

###############

sample_sheet = pd.read_csv('/scratch/henrygeo/SampleIDs.csv')
SAMPLES = sample_sheet['Sample'].tolist()

rule star_index:
    input:
      fasta="/project/6004275/genomes/Ipomoeanil_genomic.fasta",
      gtf="/project/6004275/genomes/Ipomoeanil_genomic.gtf"
    output:
      directory("Ipomoeanil_genomic")
    message:
      "Testing STAR index"
    conda: 
        "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 8000,
       time = '01:00:00'
    threads: 4
    shell:
      "STAR --runMode genomeGenerate --runThreadN 4 --genomeSAindexNbases 13 --genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"

rule star_align_pe_multi:
    input:
        fq1=expand("/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R1.fastq", sample = SAMPLES),
        fq2=expand("/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R2.fastq", sample = SAMPLES)  #optional
        idx="/scratch/henrygeo/Ipomoeanil_genomic/", # path to STAR reference genome index
    output:
        # see STAR manual for additional output files
        sam="/scratch/henrygeo/star/pe/{sample}/Aligned.out.sam",
        log="/scratch/henrygeo/star/pe/{sample}/Log.out",
    conda: 
        "envs/mapping.yaml"
    log:
        "logs/star/pe/{sample}.log",
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 28000,
       time = '01:00:00'
    threads: 8  
    script:
        "Scripts/Star_align_pe.py"