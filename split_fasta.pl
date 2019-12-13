#!/usr/bin/perl
use warnings;
use strict;

my ($file,$count) = @ARGV;

my $total = count_entries($file);

my $count_per_file = int($total/$count);

if ($total % $count) {
		++$count_per_file;
}

warn "Total seqs = $total For $count files we need $count_per_file\n";

split_file($file,$count_per_file);


sub split_file {
		my ($file,$count_per_file) = @_;

		my $prefix = $file;
		$prefix =~ s/\.fa//;

		my $seq;
		my $file_count = 0;

		open (IN,$file) or die $!;

		my $local_count = 0;

		warn "Count per file is $count_per_file\n";

		while (<IN>) {
				if (/^>/) {
						# We're at a new sequence

						if ($local_count % $count_per_file == 0) {
								my $index = ($local_count/$count_per_file)+1;
								if ($file_count) {
										close OUT or die $!;
								}
								++$file_count;

								warn "Written $local_count sequences to $index files\n";
#								sleep(1);

								open (OUT ,'>',"${prefix}_${index}.fa") or die $!;
						}
						
						print OUT $seq if ($seq);
						++$local_count;
						$seq = $_;
				}
				else {
						$seq .= $_;
				}
		}

		print OUT $seq;

		close OUT or die $!;



}


sub count_entries {
		my ($file) = @_;

		my $count = 0;

		open (IN,$file) or die $!;

		while (<IN>) {
				++$count if (/^>/);
		}

		close IN;

		return $count;
}
