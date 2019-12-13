#!/usr/bin/perl
use warnings;
use strict;

open (IN,'c_elegans.PRJNA13758.WS242.annotations.gff3') or die $!;

my %ids;

while (<IN>) {
		next if (/^#/);
	
		my ($chr,undef,$type,$start,$end,undef,$strand,undef,$annot) = split(/\t/);

		next unless ($type eq 'exon');

		next if (/Parent=(CDS|Pseudogene)/);

		if ($start > $end) {
				my $temp = $start;
				$start = $end;
				$end = $temp;
		}

		if ($annot =~ /Transcript:([\w\.]+)/) {
				my $id = uc($1);
				next if (exists $ids{$id});
				$ids{$id} = 1;
#				warn "Found ID '$id'\n";
		}
		else {
				warn "No transcript id from $annot\n";
				next;
		}
}

open (OUT,'>','rt_pcr_ids.txt') or die $!;

foreach my $id (sort keys %ids) {
		print OUT $id,"\n";
}

close OUT or die $!;
