nextflow.enable.dsl = 2

//Define parameters
params.seq_id = 'NC_003663.2'

//Build Efetch URL
process buildEfetchUrl {
    output:
    stdout

    script:
    """
    echo "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${params.seq_id}&rettype=gb&retmode=text"
    """
}
/*Create the process for receiving the GenBank file.
To download the file, we use the URL created in the previous process.
As a result, we get the GenBank file saved in the output folder.
*/
process downloadGenebank {
    publishDir "output", mode: 'copy', overwrite: true
    input:
    val (url)

    output:
    path "${params.seq_id}.gb"

    script:
    """
    wget -O ${params.seq_id}.gb "$url"
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
    buildEfetchUrl | downloadGenebank | fileSplitting
}

