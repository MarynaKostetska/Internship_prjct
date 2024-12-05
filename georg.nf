nextflow.enable.dsl = 2

params.accession = "NC_003663.2"
params.seq_id = "NC_003663.2"
params.out = "${launchDir}/output"
params.storeDir="${launchDir}/cache"


process downloadAccession {
	storeDir params.storeDir
	input:
		val accession
	output:
		path "${params.accession}.fasta1"
		
	script:
	"""
	wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${params.accession}&rettype=fasta&retmode=text" -O ${params.accession}.fasta1
	
	"""
}



process download_seq_id {
	storeDir params.storeDir
	input:
		val seq_id
	output:
		path "${params.seq_id}.gb"
		
	script:
	"""
	
	wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${params.seq_id}&rettype=gb&retmode=text" -O ${params.seq_id}.gb
	"""
}



/* Splitting a genbank file into a FASTA and GFF3 file using the genbank_to script.
For a FASTA file, use the -n option, which displays the nucleotide sequence.
For a GFF3 file, we use the --gff3 option, which generates a GFF3 file containing 9 columns
(seqid, source, type, start, end, score, strand, phase, attributes).
The 2 files obtained as a result of the process are also saved to the output folder.
*/
process fileSplitting {
    publishDir "output", mode: 'copy', overwrite: true

    input:
    path gb

    output:
	path "${params.seq_id}.fasta"
    path "${params.seq_id}.gff"
	
	script:
	"""
    genbank_to -g ${gb} -n ${params.seq_id}.fasta
    genbank_to -g ${gb} --gff3 ${params.seq_id}.gff
	"""
}



workflow {
	accession_channel = Channel.from(params.accession)
    download_accession_channel = downloadAccession(accession_channel)

	seq_id_channel = Channel.from(params.seq_id)
    download_seq_id_channel = download_seq_id(seq_id_channel)
	
	genebank_to_channel = fileSplitting(download_seq_id_channel)
	
}