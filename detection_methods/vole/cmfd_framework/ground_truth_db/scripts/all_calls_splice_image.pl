#!/usr/bin/perl -w
use strict;

# Executes `find_snippet_positions.pl` for all images in the database. 
# Effectively, this creates calls for determining the positions of the snippets
# in the copy move database. These positions have to be inserted (currently
# manually, unfortunately) in the image properties hash in db_setup.pl and
# serve as a starting point for creating custom ground-truth test cases.
# 
# note that log output is redirected to '/tmp/splice_log_$i.txt', maybe you
# want to change this

use FindBin;
push @INC, $FindBin::RealBin; # for 'do'

#my $DOWARN = 1;
#BEGIN { $SIG{'__WARN__'} = sub { warn $_[0] if $DOWARN } }
{ package SpliceConfig; do 'splice_configs.pl' }
my %splice_cfg = %SpliceConfig::splice_config;

if (scalar @ARGV < 2) {
	print "usage: $0 <output_dir> <config>\n";
	print "where <config> is one of {", join(", ", sort keys %splice_cfg), "}\n";
	exit(1);
}

my $output_dir = shift;
my $config = shift;

$| = 1;
# iterate over all images
for (my $i = 1; $i < 49; $i++) {
	my $log_file = "/tmp/splice_log_$i.txt";
	unlink $log_file if (-e $log_file);

#	print `$FindBin::RealBin/splice_image.pl $output_dir $config $i >> $log_file`;
	print "$FindBin::RealBin/splice_image.pl $output_dir $config $i >> $log_file\n";
}

