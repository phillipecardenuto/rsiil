#!/usr/bin/perl -w

use strict;
use Carp qw(cluck);

# creates ground truth using different splicing scenarios. Optionally, the
# splicing is not explicitly performed, but only the commands are printed to
# stdout.

# Dump stack trace on warnings (typically commented out)
#BEGIN { $SIG{'__WARN__'} = sub { warn $_[0]; cluck;} }

use File::Path qw(mkpath);
use File::Find;
use File::Spec::Functions;
use FindBin;
push @INC, $FindBin::RealBin; # for 'do'

my $debug = 0;

# Konfigurationsdateien einlesen
{ package Config; do 'db_setup.pl' }
{ package SpliceConfig; do 'splice_configs.pl' }

my $failure = 0;
if ((length($config::vole) == 0) || ($Config::vole eq '/vole/binary/including/full/path')) {
	print STDERR "\n\nplease edit the variable \$vole in the configuration file db_setup.pl - it should point to your vole binary, including the full program path\n";
	$failure = 1;
}
if ((length($config::vole) == 0) || ($Config::db_root eq '/path/to/the/dataset/')) {
	print STDERR "\n\nplease edit the variable \$db_root in the configuration file db_setup.pl - it should point to your dataset, sucht that appending the entries in \%files (also in db_setup.pl) points to the images.\n";
	$failure = 1;
}
exit(1) if ($failure);



if ($debug) {
	print<<EOT
	vole:           $Config::vole,
	splice_params:  $Config::splice_command_params
	splice_params2: $Config::splice_command_params2
EOT
;
}

if (scalar @ARGV < 3) {
	my $t1 = join(", ", sort keys %SpliceConfig::splice_config);
	my $t2 = scalar keys %Config::files;

	print <<EOF
usage: $0 [-c] <output_dir> <config> <numbers of images>
where
  -c is an optional argument to only create the splicing commands and print the
	 commands on stdout.  If -c is not set, the splicing commands are directly
	 executed.
  <output_dir> is the base directory for storing the spliced images
  <config> is one of {$t1}
  <numbers of images> denotes arbitrary many integer arguments, containing the
	 numbers of the images that should be spliced. Currently, numbers between 1
	 and $t2 should work.
EOF
;exit(1);
}

# check if -c is given:
my $commands_are_only_printed = 0;
@ARGV = grep {
	if ($_ =~ /^-c$/i) {
		$commands_are_only_printed = 1;
		undef;
	} else {
		$_;
	}
} @ARGV; # is the argument -c somewhere? remove it from argv, register that it has been there

my $output_dir = shift;
my $config = shift;

if (!defined $SpliceConfig::splice_config{$config}) {
	print STDERR "ERROR: config $config is not defined, should be one of {", join(", ", sort keys %SpliceConfig::splice_config), "}\n";
	exit(1);
}

my $cfg = $SpliceConfig::splice_config{$config};

# determine which params are set in this config
my @used_params = ();
my @used_command_params = ();
my @used_command_params2 = ();
for (my $i = 0; $i < scalar @Config::splice_params; $i++) {
	push @used_params, $Config::splice_params[$i] if (defined $cfg->{$Config::splice_params[$i]."_min"});
	push @used_command_params, $Config::splice_command_params[$i] if (defined $cfg->{$Config::splice_params[$i]."_min"});
	push @used_command_params2, $Config::splice_command_params2[$i] if (defined $cfg->{$Config::splice_params[$i]."_min"});
}

# double-check config for these params
foreach my $sparam(@used_params) {
	$cfg = fix_config($cfg, $sparam);
}

my $special_param_string = "";
my $ignore_all_snippets = 0;
foreach my $special_param(keys %Config::non_range_params) {
	if ((!defined $cfg->{$special_param}) or (!defined $Config::non_range_params{$special_param}{$cfg->{$special_param}})) {
		print STDERR "ERROR: parameter $special_param is missing or wrong in config $config, see db_setup.pl for admissible values\n";
		exit(1);
	}
	$special_param_string .= " " . $Config::non_range_params{$special_param}{$cfg->{$special_param}};
	if (($special_param eq 'gt_type') && ($cfg->{$special_param} eq 'none')) { # Special case: just change the input image, do not splice anything
		$ignore_all_snippets = 1;
	}
}

foreach my $bool_flag(keys %Config::bool_flags) {
	if (defined $cfg->{$bool_flag}) {
		$special_param_string .= " " . $Config::bool_flags{$bool_flag};
	}
}

# create output directory
if ((-e $output_dir) && (!-d $output_dir)) {
	print STDERR "ERROR: $output_dir exists, but is no directory, aborted\n";
	exit(1);
}

create_directory_or_print($output_dir, $commands_are_only_printed);

foreach my $number (@ARGV) {
	if ($number =~ /\D/) { print "argument $number appears to be no number, skipped.\n"; next; }
	if ($number < 1) { print "argument $number is smaller than 1, skipped.\n"; next; }
	if ($number > scalar keys %Config::files) { print "argument $number is larger than the number of known files (", scalar keys %Config::files, "), skipped.\n"; next; }

	# insert path in all relevant images
	$Config::files{$number} = expand_config_paths($Config::files{$number}, "orig");
	$Config::files{$number} = expand_config_paths($Config::files{$number}, "snippets");
	$Config::files{$number} = expand_config_paths($Config::files{$number}, "snippets_alpha");

	create_directory_or_print(catfile($output_dir, $number), $commands_are_only_printed);

	# enter loop: iterate over every range; generate & print all commands
	my $static_params =
		"--file_prefix $Config::files{$number}{prefix} " .
		"--output " . catfile($output_dir,$number) . " " .
		"--insert_positions " . $Config::files{$number}{insert_positions} . " " .
		"--source_positions " . $Config::files{$number}{source_positions} . " " .
		"--orig $Config::files{$number}{orig} " .
		"--snippet " . join(',', @{$Config::files{$number}{snippets}}) . " " .
		"--alpha_snippet " . join(',', @{$Config::files{$number}{snippets_alpha}}) . " " .
		$special_param_string
	;
	if ($ignore_all_snippets == 1) {
		$static_params =
			"--file_prefix $Config::files{$number}{prefix} " .
			"--output " . catfile($output_dir,$number) . " " .
			"--orig $Config::files{$number}{orig} " .
			$special_param_string
		;
	}

	dfs_param_list($commands_are_only_printed, $static_params, 0, "");
}

# method calls itself recursively. On the n-th recursion level, the n-th
# parameter is incremented. Terminating condition: If the recursion reaches the
# bottom (i.e. higher iter_depth >= number of parameters), the command is
# printed
sub dfs_param_list {
	my ($commands_are_only_printed, $static_params, $iter_depth, $params_accu) = @_;

	if ($iter_depth >= scalar @used_params) {
		my $command = "$Config::vole splice $static_params $params_accu";
		if ($commands_are_only_printed) {
			print "$command\n";
		} else {
			print `$command`;
		}
		return;
	}

	my $min = $used_params[$iter_depth]."_min";
	my $max = "$used_params[$iter_depth]_max";
	my $step = "$used_params[$iter_depth]_step";

	for (my $i = $cfg->{$min}; $i <= $cfg->{$max}; $i += $cfg->{$step}) {
		my $new_accu = $params_accu . " --$used_command_params[$iter_depth] --$used_command_params2[$iter_depth] $i";
		dfs_param_list($commands_are_only_printed, $static_params, $iter_depth + 1, $new_accu);
	}
}

# checks whether the configuration hash is halfways sane (does not test all
# possible cases, though), and ad-hoc fixes minor issues.
sub fix_config {
	my ($cfg, $prefix) = @_;
	my $min = "${prefix}_min";
	my $max = "${prefix}_max";
	my $step = "${prefix}_step";
	if (defined $cfg->{$min}) {
		if ( ($cfg->{$min} =~ /\D/) || ($cfg->{$max} =~ /\D/) || ($cfg->{$step} =~ /\D/)) { # configuration entries other not from 0 to 9?
			print STDERR "ERROR: Configuration options for $prefix contain non-numeric arguments; removed from the option list\n";
			undef $cfg->{$min};
		}
		if ((!defined $cfg->{$max}) || (!defined $cfg->{$step})) { # other keys not defined? do single experiment
			$cfg->{$max} = $cfg->{$min};
			$cfg->{$step} = $cfg->{$min};
		}
		if (($cfg->{$max} < $cfg->{$min}) || ($cfg->{$step} <= 0)) { # other oddities? do single experiment
			$cfg->{$step} = 1;
			$cfg->{$max} = $cfg->{$min};
		}
	}
	return $cfg;
}

# puts the path name in front of all file names that do not start with '/' or
# '../'
# (such file names are assumed to have already a path prefix contain either
# absolute or relative paths; (mental note: needs adjustment for Windows users).
sub expand_config_paths {
	my ($cfg, $key) = @_;
	my $file_names = $cfg->{$key};
	if (ref($file_names) eq "ARRAY") {
		my @new_array = map {
			my $tmp = $_;
			if (file_name_is_absolute($tmp)) {
				$tmp;
			} else {
				catfile($cfg->{path}, $tmp);
			}
		} @$file_names;

		$cfg->{$key} = \@new_array;
	} else { # file name is only a scalar - no array handling required
		if (!file_name_is_absolute($file_names)) {
			$cfg->{$key} = catfile($cfg->{path},$file_names); # no path prefix? put the image path in front
		}
	}
	return $cfg;
}

sub create_directory_or_print {
	my $dir = shift;
	my $commands_are_only_printed = shift;
	return if (-d $dir); # exit if it exists

	if ($commands_are_only_printed == 1) {
		print "mkdir $dir\n";
	} else {
		if (scalar mkpath($dir, 0) == 0) {
			print STDERR "ERROR: unable to create $dir, aborted.\n";
			exit(1);
		} else {
			print "# created directory $dir\n";
		}
	}
}
