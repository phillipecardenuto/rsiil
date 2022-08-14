#!/usr/bin/perl -w

use strict;

# separates the alpha channel out of a png file (since this is apparently not
# possible with the current version of OpenCV [1]), and stores it as a separate
# gray-level file.
# It uses imagemagick (more precisely, the imagemagick command line tool
# 'convert') to extract the channel. Note that in order to use this script,
# imagemagick must be installed on your system and the convert command must be
# in your system path. If this is not the case, consider to edit the 'convert'
# call below.
#
# Pass as argument the directory to crawl for "suitable images". A suitable
# image contains either the prefix 'snippet', or a pattern like
# '*_<number>.png' or '*_<number>_trim.png'
# The alpha channel is in a file that has the same name as the source file,
# only that '.png' is replaced by '_alpha.png'.
# Note that the alpha channel is required for all operations with the database.
# Thus, extracting the alpha channels should always be the first operation when
# working with the dataset.

# [1] Meanwhile, openCV supports alpha channels in PNG files. However, we did
#     not update the code yet.

use File::Find;

my $override = 0;

if (scalar @ARGV != 1) {
	print STDERR "usage: $0 <root_directory>\nextracts the alpha-channel from all files below <root_directory> that are named \"snippet\"\n";
}
my $root_dir = shift;

find({wanted => \&wanted, follow => 0, no_chdir => 0}, $root_dir);

sub wanted {
	# must be a file, name must start with "snippet", and must contain a ".". Then convert it.
	# alternately: endings on "_<number>.png" are also converted
	-f && ((/^snippet.*\..+$/i) || (/.*_\d.png$/i) || (/.*_\d_trim.png$/i)) && do {
		return if /_alpha\./;
		# change snippet*.* to snippet*_alpha.*
		my $file = $_;
		my $file2 = $file;
		my $path = $File::Find::name;
		my $path2 = $path;
		$file =~ s/^(.+)\.(.+?)$/$1_alpha.$2/;
		$path =~ s/^(.+)\/$file2$/$1\/$file/g;
		if ((-e $path) && (!$override)) {
			print STDERR "   found: $path2\n   target: $path\n  target file already exists, skipped - set \$override=1 to change this behavior\n";
			return;
		}
		print "found it: ", $path2, ", replaced: $file, full: $path\n";
		print "convert $path2 -channel A -negate -separate -depth 8 -type grayscale $path\n";
		`convert $path2 -channel A -negate -separate -depth 8 -type grayscale $path`;
	}
}

