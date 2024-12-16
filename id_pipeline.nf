nextflow.enable.dsl = 2

//Define parameters
params.seq_id = 'KU710884.1' //'KU710884.1'(HVB), 'MH370367.1' (Salmonella phage S114),'AH002297.2', 'NC_003663.2'
params.out = "${launchDir}/output"
 
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
    publishDir params.out
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
The 2 files obtained as a result of the process are also saved to the output folder.*/

process fileSplitting {
    publishDir params.out

    input:
    path gb

    /* 'emit' is a directive that allows to declare the output of a process with a specific name.
    This means that you can access this source data in a workflow using this name.*/
    output:
	path "${params.seq_id}.fasta", emit: fasta
    path "${params.seq_id}.gff", emit: gff
	
	script:
	"""
    genbank_to -g ${gb} -n ${params.seq_id}.fasta
    genbank_to -g ${gb} --gff3 ${params.seq_id}.gff
	"""
}

/*
process removeTitle {
    publishDir params.out
    input:
    path fasta

    output:
    path "${params.seq_id}_sequence.fasta"

    script:
    """
    awk '!/^>/ {print}' ${fasta} > "${params.seq_id}_sequence.fasta"
    """
}
*/
workflow {
    //Create URL chanenel
    efetch_url_channel = buildEfetchUrl()

    //Pass the URL to the downloadGenebank process
    genbank_channel = downloadGenebank(efetch_url_channel)

    //Pass the GenBank file in to fileSplitting process
     fasta_and_gff = fileSplitting(genbank_channel)

    //Pass only FASTA
    //no_title_fasta = removeTitle(fasta_and_gff.fasta) 

}
