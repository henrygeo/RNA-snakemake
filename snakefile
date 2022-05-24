samples = "AL1 DE11 DE13 DE14 DE15 DE20 DE28 DE31 DE35 DE36 DE39 DE7 DE9 GA1 GA3 GA4 GA5 GA7 MD1 MD12 MD13 MD21 MD26 MD27 MD30 MD34 MD40 MD45 MD5 MD53 MD54 MD55 MD56 MD57 NC1 NC15 NC21 NC22 NC24 NC27 NC29 NC35 NC42 NC8 NJ1 NJ15 NJ17 NJ28 NJ3 NJ34 NJ35 NJ41 NJ50 PA12 PA15 PA2 PA21 PA26 PA3 PA31 PA33 PA38 SC10 SC13 SC16 SC17 SC20 SC2 SC5 SC6 SC8 SC9 TN10 TN1 TN13 TN15 TN21 TN22 TN23 TN25 TN28 TN30 TN3 TN32 TN34 TN36 TN38 TN39 TN41 VA1 VA24 VA25 VA38 VA39 VA41 VA50 VA7".split()

def strbams(bams):
    bamz =[bams.format(sample=sample) for sample in samples]
    obam = " ".join([str(x) for x in bamz])
    return obam
    
rule all:
    input:  
        blam=expand("/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam", sample = samples),
        outgtf = expand("/scratch/henrygeo/Results/stringtie/{sample}.gtf", sample = samples),
        mergedgtf = "/scratch/henrygeo/Results/stringtie/merged.gtf",
        covtab = expand("/scratch/henrygeo/Results/stringtie/{sample}/e_data.ctab", sample = samples),
        
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

rule star_align_pe_first:
    input:
        fq1="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R1.fastq",
        fq2="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R2.fastq",  #optional
        idx="/scratch/henrygeo/Ipomoeanil_genomic/", # path to STAR reference genome index
    output:
        # see STAR manual for additional output files
        sam="/scratch/henrygeo/Results/star/pe/mapped/{sample}Aligned.out.sam",
        log="/scratch/henrygeo/Results/star/pe/log/{sample}Log.out",
    conda: 
        "envs/mapping.yaml"
    log:
        "logs/star/pe/{sample}.log",   
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 28000,
       time = '01:00:00'
    threads: 8
    params:
        out = "/scratch/henrygeo/Results/star/pe/{sample}",
        lgout =  "/scratch/henrygeo/Results/star/pe/{sample}"  
    shell:
       """
       ( STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.out} --outStd Log) 2> {log}
       """
rule star_align_pe_second:
    input:
        fq1=expand("/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R1.fastq", sample = samples),
        fq2=expand("/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R2.fastq", sample = samples),  #optional
        idx="/scratch/henrygeo/Ipomoeanil_genomic/", # path to STAR reference genome index
        sjdbfiles = "/scratch/henrygeo/Results/star/pe/splicejunc/{sample}SJ.out.tab",
    output:
        # see STAR manual for additional output files
        bam="/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
        log="/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}Log.out",
    conda: 
        "envs/mapping.yaml"
    log:
        "logs/star/pe/{sample}.log",   
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 28000,
       time = '01:00:00'
    threads: 8
    params:
        out = "/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}",
    shell:
       """
       ( STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.out} --sjdbFileChrStartEnd {input.sjdbfiles} \
       --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outStd Log) 2> {log}
       """
       
rule stringtie_first:
    input:
        bam = "/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
        gtf = "/project/6004275/genomes/Ipomoeanil_genomic.gtf"
    output:
        outgtf = "/scratch/henrygeo/Results/stringtie/{sample}.gtf"
    log:
        "logs/stringtie/{sample}.log"
    conda: "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 8000,
       time = '03:00:00'
    threads: 8
    shell:
        """
        (stringtie -o {output.outgtf} -v --fr -p {threads} -G {input.gtf} {input.bam}) 2> {log}
        """
rule stringtie_merge:
    input:
        gtf = "/project/6004275/genomes/Ipomoeanil_genomic.gtf",
    output:
        mergegtf = "/scratch/henrygeo/Results/stringtie/merged.gtf"
    params:
      merg = "/scratch/henrygeo/mergedsamples.txt"
    log:
        "logs/stringtie/merge/merge.log"
    conda: "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 4000,
       time = '03:00:00'
    threads: 1
    shell:
        """
        (stringtie --merge -o {output.mergegtf} -v -p {threads} -G {input.gtf} {params.merg}) 2> {log}
        """
rule stringtie_cover:
    input:
        bam = "/scratch/henrygeo/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
        gtf = "/scratch/henrygeo/Results/stringtie/merged.gtf"
    output:
        outgtf = "/scratch/henrygeo/Results/stringtie/{sample}/{sample}.gtf",
        geneab = "/scratch/henrygeo/Results/stringtie/{sample}gene_abund.tab",
        covtab = "/scratch/henrygeo/Results/stringtie/{sample}/e_data.ctab"
    params:
      opath = "/scratch/henrygeo/Results/stringtie/{sample}/"
    log:
        "logs/stringtie/cover/{sample}.log"
    conda: "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 4000,
       time = '03:00:00'
    threads: 8
    shell:
        """
        (stringtie -e -b {params.opath} -A {output.geneab} -o {output.outgtf} -v -p {threads} -G {input.gtf} {input.bam}) 2> {log}
        """