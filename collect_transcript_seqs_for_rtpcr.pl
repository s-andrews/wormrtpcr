#!/usr/bin/perl
use warnings;
use strict;

my %chromosome_seqs = read_chromosomes();

my @ids = read_ids();

my %exons = read_exons(@ids);

extract_seqs(\%chromosome_seqs,\%exons);

sub extract_seqs {
    my ($seqs,$ids) = @_;

    open (OUT,'>',"transcript_sequences.fa") or die $!;

    foreach my $gene (keys %$ids) {
				my @exons = @{$ids->{$gene}};

				unless (@exons) {
						warn "No exons found for id '$gene'\n";
						next;
				}

				my $rev = 0;
				$rev = 1 if ($exons[0]->[3] eq '-');
				
				if ($rev) {
						@exons = sort {$b->[1] <=> $a->[1]} @exons;
				}
				else {
						@exons = sort {$a->[1] <=> $b->[1]} @exons;
				}
	
				my $seq;
				my @exon_positions;

				for my $index (0..$#exons) {
						
						my $length = ($exons[$index]->[2]-$exons[$index]->[1])+1;
						if ($length > 50000) {
								die "Odd length $length in $gene";
						}
						my $exon_seq = substr($seqs->{$exons[$index]->[0]},$exons[$index]->[1]-1,$length);

						if ($rev) {
								$exon_seq = reverse($exon_seq);
								$exon_seq =~ tr/GATCgatc/CTAGctag/;
						}

						$seq .= $exon_seq;

						push @exon_positions,length($seq) unless ($index == $#exons);
				}
				
				print OUT ">$gene ",join(",",@exon_positions),"\n",$seq,"\n";

    }
		
    close OUT or die $!;
}


sub read_exons {
    my @ids = @_;

    my %exons;

    foreach my $id (@ids) {
				$exons{$id} = [];
    }

    open (IN,'c_elegans.PRJNA13758.WS242.annotations.gff3') or die $!;

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
						if (exists $exons{$id}) {
								warn "Found exon for $id from $start-$end\n";
								push @{$exons{$id}},[$chr,$start,$end,$strand];
						}
						else {
#		warn "ID '$id' not in search list\n";
						}
				}
				else {
						warn "No transcript id from $annot\n";
						next;
				}

    }

    return %exons;
}

sub read_chromosomes {

    open (IN,'c_elegans.PRJNA13758.WS242.genomic.fa') or die $!;

    my $chr;
    my %chrs;

    while (<IN>) {
				if (/>(\S+)/) {
						if ($chr) {
								warn "Length was ".length($chrs{$chr})."\n";
						}

						$chr = $1;
						warn "Reading chr $chr\n";
						next;
				}

				chomp;
				next unless ($_);
				die "Sequence but no chromosome" unless ($chr);

				$chrs{$chr} .= $_;
    }

    return %chrs;

}

sub read_ids {
    open (IN,'rt_pcr_ids.txt') or die $!;

    my @ids;

    while (<IN>) {
				chomp;
				s/[\r\n\s]//g;
				next unless ($_);

				push @ids,uc($_); 

#	warn "Found id '".uc($_)."'\n";
    }

    warn "Found a list of ".scalar @ids." ids\n";
		
    return @ids;
}
