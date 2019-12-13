#!/usr/bin/perl
use warnings;
use strict;
use LWP;

my ($seqfile) = @ARGV;
my $results_file = $seqfile;

my $temp_prefix = $seqfile;
$temp_prefix =~ s/\.fa$//;

$results_file =~ s/\.fa$/_primers.txt/;

warn "Seq=$seqfile Results=$results_file Temp=$temp_prefix\n";

my $browser = LWP::UserAgent -> new();
my $url = 'http://bilin1/cgi-bin/rtprimers.cgi';
open (RESULTS,'>',$results_file) or die $!;

read_seq_file($seqfile);

sub read_seq_file {
    
    my ($file) = @_;
    my $seq;
    my $name;
    my $splice_sites;


    open (IN,$file) or die $!;

    while (<IN>) {
				chomp;
				if (/^>/) {

						if ($seq) {
								find_primers($name,$seq,$splice_sites);
						}

						s/^>//;
						($name,$splice_sites) = split(/\s+/,$_,2);
						$seq = '';
				}
				else {
						$seq .= $_;
				}
    }
    find_primers($name,$seq,$splice_sites);

}


sub find_primers {

    my ($name,$seq,$splice_sites) = @_;

    warn "Finding primers for $name with splices $splice_sites\n";

		$splice_sites =~ s/,/ /g;

    unless ($splice_sites) {
				warn "Found no splice sites for $name - making fake ones\n";

				my @fake_sites;

				for (my $s=100;$s<length($seq);$s += 100) {
						push @fake_sites,$s
				}

				$splice_sites = join(" ",@fake_sites);
    }

    my $response = $browser->post($url,[sequence => $seq, seqname=>$name, position=>$splice_sites,gc_min=>20,gc_max=>80,tmmin=>56,tmmax=>62,tmdallow=>2,prodmin=>70,prodmax=>120,return_no=>100]);

    die "$url error: ", $response->status_line() unless $response->is_success();
    die "Weird content type at $url -- ", $response->content_type() unless $response->content_type() eq 'text/html';

    my $html =  $response->content();

		warn $html;

    if ($html =~ /No Suitable Primer Pairs found/) {
				warn "Found no primers for $name\n";
				return;
    }

    if ($html =~ /No valid splice positions found/) {
				warn "Found no primers for $name\n";
				return;
    }

    my $found = 0;

    my ($for_start,$for_end,$for_seq,$for_melt,$rev_start,$rev_end,$rev_seq,$rev_melt);

    while ($html =~ /<pre>(\d+)\s+([GATC]+)\s+(\d+)<\/pre>Tm:\s+([\d\.]+)/g) {
				if ($found == 0) {
						$for_start = $1;
						$for_seq = $2;
						$for_end = $3;
						$for_melt = $4;
				}
				elsif ($found == 1) {
						$rev_start = $1;
						$rev_seq = $2;
						$rev_end = $3;
						$rev_melt = $4;

						## Check if the primers we have are unique
						if (!(is_unique($for_seq) && is_unique($rev_seq))) {
								$found = 0;
								($for_start,$for_end,$for_seq,$for_melt,$rev_start,$rev_end,$rev_seq,$rev_melt) = (undef,undef,undef,undef,undef,undef,undef,undef);
								next;
						}
						last;
				}
				else {
						die "Found more than 2 primers\n";
				}
				++$found;
    }

    unless ($for_start && $rev_start) {
				warn "Found no primers for $name\n";
				return;
#	die $html;
    }

    print RESULTS join("\t",($name,length($seq),$for_start,$for_seq,$for_end,$for_melt,$rev_start,$rev_seq,$rev_end,$rev_melt)),"\n";
}

sub is_unique {
    my ($seq) = @_;

    warn "Checking if $seq is unique\n";

#    return (1);

    open (BLATOUT,'>',"${temp_prefix}_query.fa") or die $!;

    print BLATOUT ">query\n$seq\n";
    close BLATOUT or die $!;

    system("blat -minScore=0 -tileSize=18 -minMatch=2 -stepSize=2 -noHead -oneOff=1 c_elegans.PRJNA13758.WS242.genomic.fa ${temp_prefix}_query.fa ${temp_prefix}_blat_out.txt") == 0 or die "Failed to run blat";

    open (BLAT,"${temp_prefix}_blat_out.txt") or die $!;

    my $count = 0;
    ++$count while (my $line = <BLAT>);
		
    close (BLAT);

		unlink ("${temp_prefix}_blat_out.txt") or die $!;

    if ($count > 1) {
				warn "$seq is not unique (count = $count)\n";
				return 0;
    }

    warn "$seq is unique (count = $count)\n";
    return 1;
}
