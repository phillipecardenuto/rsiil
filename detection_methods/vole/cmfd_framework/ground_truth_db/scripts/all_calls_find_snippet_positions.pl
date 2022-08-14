#!/usr/bin/perl -w
use strict;

# Executes `find_snippet_positions.pl` for all images in the database. 
# Effectively, this creates calls for determining the positions of the snippets
# in the copy move database. These positions have to be inserted (currently
# manually, unfortunately) in the image properties hash in db_setup.pl and
# serve as a starting point for creating custom ground-truth test cases.

use FindBin;
push @INC, $FindBin::RealBin; # for 'do'

for (my $i = 1; $i < 49; $i++) {
	my @next_call = `$FindBin::RealBin/find_snippet_positions.pl $i`;

	my $log_file = "/tmp/log_$i.txt";
	unlink $log_file if (-e $log_file);
	foreach my $line(@next_call) {
		chomp($line);
		next if ($line !~ /\w/);
		print $line . " >> $log_file\n";
	}
}

