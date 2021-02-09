include: '/opt/pipelines/ngs_basic.smk'


# FIXME call_target file is needed to get chromosomes, so it must exists before running the snakemake
try:
	chroms = sorted({l.split('\t')[0] for l in open(config['call_target'], 'rt') if '\t' in l})
except FileNotFoundError:
	chroms = []
	print('Warning: target file unavailable')

wildcard_constraints:
	chrom='|'.join(chroms),
#	ref='|'.join(get_references()),

# GATK PREPROCESSING
rule gatk_MarkDuplicates:
	input:
		'{file}.bam'
	output:
		'{file}.dups.bam'
	threads: 2
	shell:
		'''gatk --java-options -Xmx4g MarkDuplicatesSpark --spark-master local[{threads}] -I {input[0]} -O {output[0]} 2> {output[0]}.log'''

rule gatk_BaseRecalibrator:
	'''Stage 1 of 2 of the base quality score recalibration (BQSR) process'''
	input:
		aln='{file}_ref={ref}.{postaln_steps}.bam',     # alignments
		ref=lambda w: get_resource(w.ref, 'sequences'), # reference genome
		kv=lambda w: get_resource(w.ref, 'dbsnp'),      # known variants
	output:
		'{file}_ref={ref}.{postaln_steps}.recal_data.table'
	shell:
		'''gatk --java-options -Xmx4g BaseRecalibrator --input {input.aln} --output {output[0]} --reference {input.ref} --known-sites {input.kv} 2> {output[0]}.log'''

rule gatk_ApplyBQSR:
	'''Stage 2 of 2 of the base quality score recalibration (BQSR) process'''
	input:
		aln='{file}_ref={ref}.{postaln_steps}.bam',              # alignments
		rec='{file}_ref={ref}.{postaln_steps}.recal_data.table', # recalibration table
		ref=lambda w: get_resource(w.ref, 'sequences'),          # reference genome
	output:
		'{file}_ref={ref}.{postaln_steps}.bqsr.bam'
	shell:
		'''gatk --java-options -Xmx4g ApplyBQSR --input {input.aln} --output {output[0]} --reference {input.ref} --bqsr-recal-file {input.rec} 2> {output[0]}.log'''


# GATK GERMLINE
rule gatk_HaplotypeCaller:
	# single sample germline call with gatk HaplotypeCaller
	input:
		aln='{file}_ref={ref}.{postaln_steps}.bam',        # sorted alignments
		alnidx='{file}_ref={ref}.{postaln_steps}.bam.bai', # alignment index
		target=config['call_target'],                      # target bed
		ref=lambda w: get_resource(w.ref, 'sequences'),    # reference genome
	output:
		'{file}_ref={ref}.{postaln_steps}.gatk_hc.g.vcf.gz'
	threads: 1
	resources:
		mem_mb=4000
	shell:
		'''gatk --java-options -Xmx4g HaplotypeCaller -R {input.ref} -I {input.aln} --intervals {input.target} -O {output[0]} -ERC GVCF --native-pair-hmm-threads {threads} 2> {output[0]}.log'''

def get_sample_list(path, suffix):
	import pandas
	ss = pandas.read_csv(path, sep='\t', index_col='SAMPLE')
	return (ss['PATH'] + suffix).to_dict()
	
rule gatk_sample_sheet:
	input:
		'{sample_sheet}.tsv',
		unpack(lambda w: get_sample_list(w.sample_sheet + '.tsv', '.' + w.params + '.g.vcf.gz')),
	output:
		'{sample_sheet,[^./]+}.{params}.gatk_gendb/sample_map.tsv'
	run:
		with open(output[0], 'wt') as f:
			for sample, path in input.items():
				print(sample, path, sep='\t', file=f)

rule gatk_GenomicsDBImport:
	# collect all single sample gvcf into a database before the joint calling
	input:
		'{file}/sample_map.tsv' #FIXME add depend on vcf files
	output:
		directory('{file}/{chrom}')
	shell:
		'''gatk --java-options -Xmx4g GenomicsDBImport --sample-name-map {input[0]} --genomicsdb-workspace-path {output} -L {wildcards.chrom} 2> {output}.log'''

rule gatk_GenotypeGVCFs:
	# joint calling 
	input:
		lambda w: get_resource(w.ref, 'sequences'),
		'{part1}_ref={ref}.{part2}/{chrom}',
	output:
		'{part1}_ref={ref}.{part2}/{chrom}.gatk_ggvcf.vcf.gz'
	threads: 1
	shell:
		'gatk --java-options -Xmx4g GenotypeGVCFs -R {input[0]} -V gendb://{input[1]} -O {output} 2> {output}.log'

rule bcftools_concat_chroms:
	input:
		expand('{{file}}/{chrom}.gatk_ggvcf.vcf.gz', chrom=chroms)
	output:
		'{file}.allchr.gatk_ggvcf.vcf.gz'
	shell:
		'bcftools concat -Oz -o {output} {input} 2> {output}.log'

