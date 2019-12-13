# wormrtpcr
The scripts used for Laetitia Chauve's paper on high throughput worm RTPCR

The purpose and order of use of the scripts is as below:

###get_all_ids.pl
This takes in a GFF annotation file (we used c_elegans.PRJNA13758.WS242.annotations.gff3 from Ensembl), and produces a non-redundant list of transcript ids found within it.

###collect_transcript_seqs_for_rtpcr.pl
This takes the set of ids from the first script along with the genome sequence for worm (c_elegans.PRJNA13758.WS242.genomic.fa) and generates a multi-fasta file of the extracted transcript sequences along with the positions of the splice junctions within them - information which will be used by the design process later

###split_fasta.pl
This is an optional script which can split the output of the collect_transcript_seqs_for_rtpcr.pl into multiple chunks to allow for parallel processing.

###design_primers.pl
This reads the transcript sequence file from collect_transcript_seqs_for_rtpcr.pl and feeds the individual sequences to the rtprimers.cgi script to generate the list of hits.  It then takes the best hit and puts this into a common output file.

###rtprimers.cgi
This is the script which does the primer design.  It is currently a CGI script which needs to be hosted on a web server, although it could be easily adapted to run as a stand-alone script.   It uses the transcript sequences to design primers which span at least one of the exon junctions (if present), filtering potential primers on GC content, melt temperature, melt tempreature differnece and product size.

The default settings for the filters are:

* Min GC = 20%
* Max GC = 80%
* Min Tm = 58
* Max Tm = 60
* Tm Diff = 0.5
* Min Product Size = 70
* Max Product Size = 120


