samples = "AL1 DE11 DE13 DE14 DE15 DE20 DE28 DE31 DE35 DE36 DE39 DE7 DE9 GA1 GA3 GA4 GA5 GA7 MD1 MD12 MD13 MD21 MD26 MD27 MD30 MD34 MD40 MD45 MD5 MD53 MD54 MD55 MD56 MD57 NC1 NC15 NC21 NC22 NC24 NC27 NC29 NC35 NC42 NC8 NJ1 NJ15 NJ17 NJ28 NJ3 NJ34 NJ35 NJ41 NJ50 PA12 PA15 PA2 PA21 PA26 PA3 PA31 PA33 PA38 SC10 SC13 SC16 SC17 SC20 SC2 SC5 SC6 SC8 SC9 TN10 TN1 TN13 TN15 TN21 TN22 TN23 TN25 TN28 TN30 TN3 TN32 TN34 TN36 TN38 TN39 TN41 VA1 VA24 VA25 VA38 VA39 VA41 VA50 VA7".split()

rule all:
    input:  
        #mergedgtf = "/scratch/henrygeo/RNAproject/Results/stringtie/merged.gtf",
        #covtab = expand("/scratch/henrygeo/RNAproject/Results/stringtie/{sample}/e_data.ctab", sample = samples),
        sf = expand("/scratch/henrygeo/RNAproject/Results/salmon_quant/{sample}/", sample = samples),
        #outgtf = expand("/scratch/henrygeo/RNAproject/Results/stringtie/first/{sample}.gtf", sample = samples),
        #stats = "/scratch/henrygeo/RNAproject/Results/cufflinks/stringtieresults.stats",
        #bam = expand("/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",sample = samples),
        #sj = "/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}SJ.out.tab", 
        #sam = expand("/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}Aligned.out.sam",sample = samples),
        #staridx = "/scratch/henrygeo/RNAproject/Results/Ipomoeanil_genomic"
        
rule star_index:
    input:
      fasta="/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.fna",
      gtf="/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
    output:
      directory("/scratch/henrygeo/RNAproject/Results/Ipomoeanil_genomic")
    message:
      "Testing STAR index"
    conda: 
        "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 8000,
       time = '01:00:00'
    threads: 8
    shell:
      "STAR --runMode genomeGenerate --runThreadN 4 --genomeSAindexNbases 13 --genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"

#rule star_align_pe_first:
#    input:
#        fq1="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R1.fastq", 
#        fq2="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R2.fastq",   #optional
#        idx="/scratch/henrygeo/RNAproject/Results/Ipomoeanil_genomic", # path to STAR reference genome index
#    output:
#        # see STAR manual for additional output files
#        sam="/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}Aligned.out.sam",
#        log="/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}Log.out",
#        sjdb = "/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}SJ.out.tab",
#    conda: 
#        "envs/mapping.yaml"
#    resources:
#       mem_mb = lambda wildcards, attempt: attempt * 28000,
#       time = '01:00:00'
#    threads: 8
#    params:
#        out = "/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}",
#        lgout =  "/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}"  
#    shell:
#       "STAR --runThreadN {threads} --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.out} --outStd Log"

#rule star_align_pe_second:
#    input:
#        fq1="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R1.fastq",
#        fq2="/project/6004275/sharing/IpomoeahedSeq2022/{sample}_R2.fastq", #optional
#        idx="/scratch/henrygeo/RNAproject/Results/Ipomoeanil_genomic/", # path to STAR reference genome index
#        sjdbfiles = "/scratch/henrygeo/RNAproject/Results/star/pe/mapped/{sample}SJ.out.tab",
#    output:
#        # see STAR manual for additional output files
#        bam="/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
#        tbam = "/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.toTranscriptome.out.bam"
#    conda: 
#        "envs/mapping.yaml"
#    resources:
#       mem_mb = lambda wildcards, attempt: attempt * 28000,
#       time = '01:00:00'
#    threads: 8
#    params:
#        out = "/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}",
#    shell:
#       "STAR --runThreadN {threads} --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.out} --sjdbFileChrStartEnd {input.sjdbfiles} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outStd Log"

rule stringtie_first:
    input:
        bam = "/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
        gtf = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
    output:
        outgtf = "/scratch/henrygeo/RNAproject/Results/stringtie/first/{sample}.gtf"
    conda: 
        "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 8000,
       time = '03:00:00'
    threads: 8
    shell:
        "stringtie -m 150 --conservative -o {output.outgtf} -v --fr -p {threads} -G {input.gtf} {input.bam}"
        
rule stringtie_merge:
    input:
        #gtf = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
        check = expand("/scratch/henrygeo/RNAproject/Results/stringtie/first/{sample}.gtf", sample = samples)
    output:
        mergegtf = "/scratch/henrygeo/RNAproject/Results/stringtie/mergednogen.gtf"
    params:
      merg = "/scratch/henrygeo/RNAproject/mergedsamples.txt"
    conda: 
        "envs/mapping.yaml"
    resources:
       mem_mb = lambda wildcards, attempt: attempt * 4000,
       time = '03:00:00'
    threads: 1
    shell:
        "stringtie --merge -c 1 -f 0.05 -o {output.mergegtf} -v -p {threads} {params.merg}"

#rule cuffmerge:
#    input:
#        refgtf = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
#        check = expand("/scratch/henrygeo/RNAproject/Results/stringtie/first/{sample}.gtf", sample = samples)
#    output:
#        merged = "/scratch/henrygeo/RNAproject/Results/cufflinks/merged.gtf"
#    params:
#        merg = "/scratch/henrygeo/RNAproject/mergedsamples.txt",
#    conda:
#        "envs/mapping.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '03:00:00'
#    threads: 4
#    shell:
#        "cuffmerge -p {threads} -g {input.refgtf} {params.merg}"
#        cuffmerge -p 4 --no-update-check TRUE -g /scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf /scratch/henrygeo/RNAproject/mergedsamples.txt 

#rule gffcompare:
#    input:
#        ref = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
#        stringtie = "/scratch/henrygeo/RNAproject/Results/stringtie/merged.gtf",
#        cufflink = "/scratch/henrygeo/RNAproject/Results/cufflinks/merged.gtf"
#    output:
#        stats = "/scratch/henrygeo/RNAproject/Results/gffcompare/stringtiecufflinks.stats",
#        mergedgtf = "/scratch/henrygeo/RNAproject/Results/gffcompare/stringtiecufflinks.combined.gtf",
#    params:
#        prefix = "/scratch/henrygeo/RNAproject/Results/cufflinks/stringtieresults"
#    conda:
#        "envs/mapping.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 4000,
#        time = '03:00:00'
#    threads: 8
#    shell:
#        "gffcompare -r {input.ref} -R -o {params.prefix} {input.stringtie} {input.cufflink}"

rule gffread:
      input:
          mer = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.gtf",
          ref = "/scratch/henrygeo/RNAproject/genomes/GCF_001879475.1_Asagao_1.1_genomic.fna"
      output:
          trfa = "/scratch/henrygeo/RNAproject/Results/gffread/assembly_transcripts.fa"
      container:
          "docker://bschiffthaler/gffread"
      resources:
          mem_mb = lambda wildcards, attempt: attempt * 4000,
          time = '00:30:00'
      threads: 1
      shell:
          "gffread -w {output.trfa} -g {input.ref} {input.mer}"

rule Salmon:
      input:
          trfa = "/scratch/henrygeo/RNAproject/asagao.fa",
          tbam = "/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.toTranscriptome.out.bam"
      output:
          sf = directory("/scratch/henrygeo/RNAproject/Results/salmon_quant/{sample}/")
      params:
          pf = "/scratch/henrygeo/RNAproject/Results/salmon_quant/{sample}"
      container:
          "docker://combinelab/salmon"
      resources:
          mem_mb = lambda wildcards, attempt: attempt * 4000,
          time = '01:00:00'
      threads: 8
      shell:
          "salmon quant -p {threads} -t {input.trfa} -l ISF -a {input.tbam} -o {params.pf}"
        
#rule stringtie_cover:
#    input:
#        bam = "/scratch/henrygeo/RNAproject/Results/star/pe/2pass/mapped/{sample}Aligned.sortedByCoord.out.bam",
#        gtf = "/scratch/henrygeo/RNAproject/Results/stringtie/merged.gtf"
#    output:
#        outgtf = "/scratch/henrygeo/RNAproject/Results/stringtie/{sample}/{sample}.gtf",
#        geneab = "/scratch/henrygeo/RNAproject/Results/stringtie/{sample}gene_abund.tab",
#        covtab = "/scratch/henrygeo/RNAproject/Results/stringtie/{sample}/e_data.ctab"
#    params:
#      opath = "/scratch/henrygeo/RNAproject/Results/stringtie/{sample}/"
#    conda: "envs/mapping.yaml"
#    resources:
#       mem_mb = lambda wildcards, attempt: attempt * 4000,
#       time = '03:00:00'
#    threads: 8
#    shell:
#        "stringtie -e -b {params.opath} -A {output.geneab} -o {output.outgtf} -v -p {threads} -G {input.gtf} {input.bam}"