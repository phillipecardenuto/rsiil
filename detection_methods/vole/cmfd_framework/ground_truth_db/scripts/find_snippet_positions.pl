#!/usr/bin/perl -w

use strict;

# Creates calls for determining the positions of the snippets for one or
# multiple images in the copy move database. These positions have to be
# inserted (currently manually, unfortunately) in the image properties hash in
# db_setup.pl and serve as a starting point for creating custom ground-truth
# test cases.

use File::Find;
use FindBin;
push @INC, $FindBin::RealBin; # for 'do'

my $debug = 0;

# Konfigurationsdateien einlesen
{ package Config; do 'db_setup.pl' }

my $failure = 0;
if ($Config::vole eq '/vole/binary/including/full/path') {
	print STDERR "\n\nplease edit the variable \$vole in the configuration file db_setup.pl - it should point to your vole binary, including the full program path\n";
	$failure = 1;
}
if ($Config::db_root eq '/path/to/the/dataset/') {
	print STDERR "\n\nplease edit the variable \$db_root in the configuration file db_setup.pl - it should point to your dataset, sucht that appending the entries in \%files (also in db_setup.pl) points to the images.\n";
	$failure = 1;
}
exit(1) if ($failure);

if ($debug) {
	my $nfiles = scalar keys %Config::files;
	print<<EOT

	vole:  $Config::vole,
	files: $nfiles images found
EOT
;
}

if (scalar @ARGV == 0) {
	print "usage: $0 <numbers of images>\n";
	exit(1);
}

foreach my $number (@ARGV) {
	if ($number =~ /\D/) { print "argument $number appears to be no number, skipped.\n"; next; }
	if ($number < 1) { print "argument $number is smaller than 1, skipped.\n"; next; }
	if ($number > scalar keys %Config::files) { print "argument $number is larger than the number of known files (", scalar keys %Config::files, "), skipped.\n"; next; }

	my %entry = %{$Config::files{$number}};
	my $num_snippets = scalar @{$entry{snippets}};
	print "\n";
	for (my $i = 0; $i < $num_snippets; $i++) {
		my $orig = "--orig $entry{path}$entry{orig}";
		my $copy = "--copy $entry{path}$entry{copy}";
		my $snippet = "--snippet $entry{path}$entry{snippets}[$i]";
		my $snippet_alpha = "--alpha_snippet $entry{path}$entry{snippets_alpha}[$i]";
		my $ground_truth = "--ground_truth $entry{path}$entry{gt_snippets}[$i]";
		my $tolerance = "";
		if (defined $entry{l2dist}) {
			$tolerance = " --l2dist";
		} else {
			$tolerance = " --l1tol $entry{l1tol}[$i]" if (defined $entry{l1tol});
		}
		my $prune = "";
		$prune = " --prune $entry{prune}" if (defined $entry{prune});
		my $opaque = "";
		$opaque = " --opaque $entry{opaque}" if (defined $entry{opaque});
		my $full_snippets_alpha = "";
		if (defined $entry{full_snippets_alpha} and $entry{full_snippets_alpha}[$i]) {
			$full_snippets_alpha = " --full_snippet_alpha $entry{path}$entry{full_snippets_alpha}[$i]";
		}
		print "$Config::vole gt $orig $copy $snippet $snippet_alpha $ground_truth$full_snippets_alpha$tolerance$prune$opaque\n";
	}
}

print "\n";
