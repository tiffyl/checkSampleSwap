#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// HELP MESSAGE
if ( params.help ) {
    help = """sampleswaptest.nf:
             |
             |Required arguments:
             |  --nanovcf           Path to nanopore VCF.
             |  --strandseqdir      Path to directory with Strand-Seq bam files.
             |
             |Optional arugments:
             |  --chr <string>      Chromosome of interest. [default: ${params.chr}]
             |  --window <int>      Size of window to break chromosome. [default: ${params.window}]
             |  --keepbinvcf        Keep binned vcf files. [default: ${params.keepbinvcf}]
             |  --normalize <int>   Number of SNPs to compare. [default: all available]
             |  --binvcfdir <path>  Path to existing binned vcf files (will skip binning in pipeline and use existing files instead).
             |
             |Optional Nextflow arugments:
             |  -N <str>            Email to receive workflow notification.
             |  -resume             Resume nextflow execution from previous execution.
             | 
             """.stripMargin()

    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// REQUIRED
if ( ! params.nanovcf ) {
    error("Error: Please specify the if data is paired or single ended the --paired parameter [true/false].")
}

if ( ! params.strandseqdir ) {
    error("Error: Please specify the if data is paired or single ended the --paired parameter [true/false].")
}

// PROCESS
process bin_vcf {
    // bin VCF into bin size
    container "${projectDir}/singularity/bcftools.sif"

    input:
        val(windowsize)
        path(nanovcf)

    output:
        path("chr*.vcf.gz*"), optional:true, emit: subvcfs
        stdout

    script:
    """
    chrlength=\$(cat /projects/lansdorp/references/hg38/hg38.chrom.sizes | head -24 | grep -w ${params.chr} | cut -f2)
    readarray -t regions < <(bash ${projectDir}/scripts/bin_genomes.sh ${windowsize} ${params.chr} 1 \$chrlength)
    index=0

    while [[ \$index -lt \${#regions[@]} ]]; do
        region=\${regions[\$index]}

        bcftools view -Oz -r \$region -o \$region.vcf.gz ${nanovcf[0]}
        bcftools index \$region.vcf.gz

        snps=\$(bcftools +counts \$region.vcf.gz | grep "Number of SNPs:" | grep -o "[0-9]\\+")
        if [[ \$snps -lt 1 ]]; then
            rm \$region.vcf.gz*
        elif [[ \$snps -gt ${params.maxsnps} ]]; then
            rm \$region.vcf.gz*

            read -r chr start end < <(echo \$region | sed 's/\\(:\\|-\\)/ /g')

            times=\$(( ( \$snps / ( ${params.maxsnps} * 2 ) + 1) * 2 ))
            window=\$(( ( \$end - \$start ) / \$times + 1 ))

            readarray -t newregions < <(bash ${projectDir}/scripts/bin_genomes.sh \$window \$chr \$start \$end)
            regions=(\${regions[@]} \${newregions[@]})
        fi
        index=\$((index + 1))
    done 
    """
}

process out_binvcf {
    input:
        path(binnedvcfs)

    script:
    """
    mkdir -p ${launchDir}/bin_vcf/
    rsync -arL --include="*vcf.gz*" ./ ${launchDir}/bin_vcf/
    """

}

process mpileup_strandseq {
    container "${projectDir}/singularity/bcftools.sif"

    memory "2GB"

    input:
        tuple val(region), path(subnanovcf)
        path(strandseqdir)

    output:
        path("*mpileup.vcf.gz"), emit: mpileups
        path("*mpileup.vcf.gz.csi"), emit: mpileupsindex

    script:
    """   
    bcftools mpileup --no-reference -Oz -Q 30 -R ${subnanovcf[0]} -o sseq_${region}.mpileup.vcf.gz ${strandseqdir}/*bam
    bcftools index sseq_${region}.mpileup.vcf.gz
    """
}

process concat_mpileup {
    container "${projectDir}/singularity/bcftools.sif"

    input:
        val(mpileupList)
        val(mpileupindexList)

    output:
        path("allmpileup.vcf.gz"), emit: allmpileup

    script:
    """
    bcftools concat ${mpileupList.join(" ")} | bcftools view -Oz -v snps -o allmpileup.vcf.gz
    """
}

process discordance {
    container "${projectDir}/singularity/py310_viz.sif"

    input:
        path(allmpileup)
        path(nanovcf)
    
    output:

    script:
    """
    discordance=\$(python ${projectDir}/scripts/comparegt.py ${nanovcf[0]} ${allmpileup} ${params.chr} ${params.normalize})

    echo ${params.nanovcf} ${params.strandseqdir} \$discordance >> ${launchDir}/discordance.txt
    """
}


// WORKFLOW
workflow {
    Channel.fromPath("${params.nanovcf}*")
        .collect()
        .set{ nanovcf_ch }

    Channel.fromPath("${params.strandseqdir}")
        .collect()
        .set{ strandseqdir_ch }

    if ( ! params.binvcfdir ) {
        bin_vcf(params.window, nanovcf_ch)
        
        if ( params.keepbinvcf ) {
            out_binvcf(bin_vcf.out.subvcfs.collect())
        }

        bin_vcf.out.subvcfs
            .flatMap{ it }
            .map{ file -> [ file.baseName.split("\\.")[0] , file ]}
            .groupTuple()
            .set{ binnedvcf_ch }

        mpileup_strandseq(binnedvcf_ch, strandseqdir_ch)
    } 
    else {
        Channel.fromPath("${params.binvcfdir}/*.vcf.gz*")
            .map{ file -> [ file.baseName.split("\\.")[0], file ] }
            .groupTuple()
            .set { binnedvcf_ch }

        mpileup_strandseq(binnedvcf_ch, strandseqdir_ch)
    }
    
    concat_mpileup(mpileup_strandseq.out.mpileups.collect(), mpileup_strandseq.out.mpileupsindex.collect())
    discordance(concat_mpileup.out.allmpileup, nanovcf_ch)
}
