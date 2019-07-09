#!/usr/bin/perl

=head1 NAME

seqdoc.pl - Perl cgi script to align and subtractively compare two sequence
chromatograms in ABI (Applied Biosystems) format.

=head1 SYNOPSIS

seqdoc.pl?ref_seq=<reference sequence file>&test_seq=<test_sequence file>

Designed to be called as a cgi script from a web form with the two filenames
supplied as parameters.

=head1 DESCRIPTION

This script does the following in order:

=over 3

=item *

Reads in the two sequence files

=item *

Normalizes the peak heights for each trace

=item *

Aligns the two traces

=item *

Calculates a difference profile between the traces

=item *

Processes the difference profile to highlight differences characteristic of
single bases

=item *

Outputs aligned images of the two input sequences and the difference profile
in html format

=back
	
All processing is performed automatically, with no user-defined parameters other
than the sequence file names. Settings can be modified by altering values
defined in the code if desired.

The program requires write access to a temporary directory (by default /var/www/html/Temp/) to store the image files between their generation and display on the web page.

=head1 WEBSITE

http://research.imb.uq.edu.au/seqdoc/

=head1 AUTHOR

Written by Mark Crowe - marklcrowe@gmail.com

The Institute for Molecular Bioscience, The University of Queensland, Brisbane, Queensland, Australia.

=head1 COPYRIGHT

Copyright 2005, Mark Crowe.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (http://www.gnu.org/licenses/gpl.html)
for more details.

=cut




use strict;
#use warnings;
use CGI;
use ABI;
use GD::Graph::lines;

my $cgi = new CGI;

# Create a list of the four channels to save on having to keep retyping them
# in loops
our @trace_list = qw( a_trace g_trace c_trace t_trace );
our $temp_dir = '/var/www/html/Temp/';

# Get input data parameters from the web page
my $ref_seq = $cgi->param('ref_seq');
my $test_seq = $cgi->param('test_seq');

# Create temporary files to store the uploaded chromatograms
### Might need to modify to generate 'unique' names each time ###
# Reference sequence
open REF, ">".$temp_dir.$$."ref_seq";
my $ref_fh = $cgi->upload('ref_seq');
print REF while <$ref_fh>;
# Need to close files, otherwise ABI.pm can't access them properly
close REF;
# And test sequence
open TEST, ">".$temp_dir.$$."test_seq";
my $test_fh = $cgi->upload('test_seq');
print TEST while <$test_fh>;
# Need to close files, otherwise ABI.pm can't access them properly
close TEST;
# Mac uploads put in extra header info. Delete it, since all it seems to do
# is break things
demac();

# Extract the data from the chromatogram files. Get it as a hash reference
my $ref_data = read_chromatogram($temp_dir.$$.'ref_seq');
my $test_data = read_chromatogram($temp_dir.$$.'test_seq');
# Extract the sequence names
my $ref_name = $ref_data->{name};
my $test_name = $test_data->{name};

# Normalize peak heights to allow valid comparisons between runs
# By passing a hash reference, the changes made in the subroutine will act on
# the original hash (which is what we want in this case)
normalize($ref_data);
normalize($test_data);

# Now align test trace to reference trace. Reference trace doesn't change.
get_best_align($ref_data, $test_data);

# Calculate differences between traces
my ($align_length, $diffs) = differences($ref_data, $test_data);

# Finally generate output images - need to create temporary files, as can't
# include image data directly in html
# Pass details about image size and scale to the subroutine, also output file
my $ref_image = get_image($ref_data, 0, 2000, 200, $temp_dir.$$."ref_image.png", $align_length, 1);
my $test_image = get_image($test_data, 0, 2000, 200, $temp_dir.$$."test_image.png", $align_length, 1 );
my $diff_image = get_image($diffs, -500, 500, 300, $temp_dir.$$."diff_image.png", $align_length, 0);


# Print out page
print $cgi->header;
print "<head><title>Chromatogram comparison</title></head>";
print "<body>";
print "Reference sequence";
if ($ref_data->{name}) {print " - ".$ref_data->{name}."<br>"}
else {print " - ".$ref_seq."<br>"}
print "<img src=\"/Temp/".$$."ref_image.png\"><br>\n";
print "<img src=\"/Temp/".$$."diff_image.png\"><br>\n";
print "Test sequence";
if ($test_data->{name}) {print " - ".$test_data->{name}."<br>"}
else {print " - ".$test_seq."<br>"}
print "<img src=\"/Temp/".$$."test_image.png\"><br>\n";
print "</body></html>";

# Leave output files - they are deleted automatically each night


sub error_end {
	# Deals with errors in input file format
	my $error = shift;
	my %error_reply = (abi_format  => 'One of the files does not seem to be a correctly formatted ABI file.',
	                   trace_error => 'Unable to extract trace data from one of the chromatograms.',
										 basecalls   => 'Unable to extract basecall positions from one of the input files.',
										 sequence    => 'Unable to extract called sequence from one of the input files.',
										 short_seq   => 'One of the sequences is too short for a valid comparison to be made',
										);
	my $reply = $error_reply{$error};
	# Print out page
	print $cgi->header;
	print "<head><title>Chromatogram comparison error</title></head>";
	print "<body>";
	print "<h1>Error report</h1>\nAn error has occurred. The error is :-  <br><pre>$reply</pre>\nPlease check that both your input files are valid ABI chromatograms.";
	print "</body></html>";
	exit;
}

										 

sub read_chromatogram {
	# Reads in data from the chromatogram files and returns it as a hash reference
	# Provide filename as the argument to the function call
	my $file = shift;
	# Create a new ABI object pointing to that file
	my $abi = ABI->new(-file=>$file) or error_end('abi_format');
	# Get the sample name from the file
	my $sample_name = $abi->get_sample_name();
	#Now get the four traces
	my @trace_a = $abi->get_trace("A") or error_end('trace_error');
	my @trace_g = $abi->get_trace("G") or error_end('trace_error');
	my @trace_c = $abi->get_trace("C") or error_end('trace_error');
	my @trace_t = $abi->get_trace("T") or error_end('trace_error');
	# Remove blank sequence from end of trace
	# Blank sequence often has noise, so can't use equal to zero. Instead have
	# threshold cutoff. Tried with various values, even as high as 500 seems to
	# work and not delete real sequence, while as low as 5 still removes most
	# blank sequence. 50 therefore seems to be a good compromise.
	while ($trace_a[-1] < 50 && $trace_g[-1] < 50 && $trace_c[-1] < 50 && $trace_t[-1] < 50) {
	pop(@trace_a); pop(@trace_g); pop(@trace_c); pop(@trace_t)}
	
	my @base_pos = $abi->get_base_calls() or error_end('basecalls');
	my $sequence = $abi->get_sequence() or error_end('sequence');
	# Length is number of datapoints in the sequence
	my $seq_length = $#trace_a;
	# Put all the data together into a hash to return it
	my %return_hash = (name       => $sample_name,
										 a_trace    => \@trace_a,
										 g_trace    => \@trace_g,
										 c_trace    => \@trace_c,
										 t_trace    => \@trace_t,
										 seq_length => $seq_length,
										 sequence   => $sequence,
										 base_pos   => \@base_pos,
										);
	return \%return_hash;
}

sub normalize {
	# Takes the raw data from the chromatograms and normalizes each value to a 
	# mean of 100 with respect to all values within 500 datapoints either side
	# N.B. There is a special case for points within 500bp of the start.
		
	my $trace_data = shift;
	# Get the length of good trace data
	my $trace_length = $trace_data->{seq_length};
	# Ensure that the sequence is long enough for normalization to not give any 
	# errors. A shorter sequence probably won't be a lot of use anyway. Use
	# 1100 as first 500 + last 500 plus a short middle buffer. Actual value has
	# no meaning beyond preventing errors
	error_end('short_seq') unless $trace_length > 1100; 
	
	# Need to create new reference arrays, since normalization alters the existing
	# values in the hash
	my %orig_values;
	for (@trace_list) {
		for my $index (0..$trace_length) {
			$orig_values{$_}[$index] = $trace_data->{$_}[$index];
		}
	}
	
	# Do the middle of the trace first - from #500 to 500 before the end
	my $totalsum = 0;
	# Set the normalization sum before starting
	for my $index (0..999) {
		for (@trace_list) {
			$totalsum += $orig_values{$_}[$index];
		}
	}
	for my $datapoint (500..($trace_length - 501)) {
		# Normalize to 100. Divide by sum of all values, multiply by number of 	
		# values, and multiply by 100;
		# First add values for datapoint +500
		for (@trace_list) {
			$totalsum += $orig_values{$_}[$datapoint + 500];
		}
		for (@trace_list) {
			# Blank sequence can cause problems through division by zero errors.
			# Deleting trailing blank sequence helps, but just put in a default value
			# for totalsum in case of problems.
			$totalsum = 1000 unless $totalsum;
			# Calculate normalized data
			$trace_data->{$_}[$datapoint] = int(($orig_values{$_}[$datapoint] / $totalsum) * 4000 * 100);
		}
		# Finally subtract values for datapoint -500
		for (@trace_list) {
			$totalsum -= $orig_values{$_}[$datapoint - 500];
		}
	}
	# Now do first 500 - special case, since can't do 500 before. Instead just 
	# take all points before. Not so critical anyway, since data quality is poor
	# at start so any mismatches will be unreliable anyway.
	# Start by initialising totalsum
	$totalsum = 0;
	for my $index (0..500) {
		for (@trace_list) {
			$totalsum += $orig_values{$_}[$index];
		}
	}
		# Blank sequence can cause problems through division by zero errors.
		# Deleting trailing blank sequence helps, but just put in a default value
		# for totalsum in case of problems.
		$totalsum = 1000 unless $totalsum;
		
	for my $datapoint (0..499) {
		# Can do 500 after though
		my $end = $datapoint + 500;
		# Normalize to 100. Divide by sum of all values, multiply by number of
		# values, and multiply by 100;
		for (@trace_list) {
			# Blank sequence can cause problems through division by zero errors.
			# Deleting trailing blank sequence helps, but just put in a default value
			# for totalsum in case of problems.
			$totalsum = 1000 unless $totalsum;
			# Calculate normalized data
			$trace_data->{$_}[$datapoint] = int(($orig_values{$_}[$datapoint] / $totalsum) * $end * 4 * 100);
		}
		# Add next value to totalsum, to keep 500 values after
		for (@trace_list) {
			$totalsum += $orig_values{$_}[$end + 1];
		}
	}
	
	# Finally the last 500 - again a special case, as can't do 500 after. Instead 
	# just take all points after. Not so critical anyway, since data quality is
	# poor at end so any mismatches will be unreliable anyway.
	# Again start by initialising totalsum
	$totalsum = 0;  
	for my $index (($trace_length-1000)..($trace_length-1)) {
		for (@trace_list) {
			$totalsum += abs($orig_values{$_}[$index]);
		}
	}
	# Identify start and finish points
	my ($first, $last) = (($trace_length-500),($trace_length-1));
	for my $datapoint ($first..$last) {
		my $start = $datapoint - 500;
		# Normalize to 100. Divide by sum of all values, multiply by number of
		# values, and multiply by 100;
		for (@trace_list) {
			# Blank sequence can cause problems through division by zero errors.
			# Deleting trailing blank sequence helps, but just put in a default value
			# for totalsum in case of problems.
			$totalsum = 1000 unless $totalsum;
			# Calculate normalized data
			$trace_data->{$_}[$datapoint] = int(($orig_values{$_}[$datapoint] / $totalsum) * ($last-$start) * 4 * 100);
		}
		# Subtract first value from totalsum, to keep to 500 point before test
		for (@trace_list) {
			$totalsum -= abs($orig_values{$_}[$start]);
		}	
	}

	return 1;
}

		
sub get_best_align {
# This does an alignment of the first 1000 datapoints using a range of offsets
# from -200 to 200. The best alignment is picked on the basis of having the
# lowest score from datapoint 200 to 1000, and is used to allow for any
# variation in start position of the two sequences.
	my ($ref, $test) = @_;
	my %scores;
	for (my $offset=-200;$offset<=200;$offset+=20) {
		# Create temporary hashes, since a hash reference is passed to the function
		# and otherwise the real hash will be modified
		my (%temp_ref, %temp_test);
		# Also need to extract hash values, since they too are references and 
		# can't be copied for the same reason.
		for (@trace_list) {
			my @temp_ref_trace = @{$ref->{$_}};
			$temp_ref{$_} = \@temp_ref_trace;
			my @temp_test_trace = @{$test->{$_}};
			$temp_test{$_} = \@temp_test_trace;
		}
		# Do a partial alignment (first 1000 datapoints)
		align (\%temp_ref, \%temp_test, $offset, 1000);
		# Work out the score for that alignment
		$scores{$offset} = get_score(200, 1000, 0, \%temp_ref, \%temp_test);
	}
	# Sort the scores to find out the lowest, and record the value of that offset
	my @offset = sort { $scores{$a} <=> $scores{$b} } keys %scores;
	# Once the best alignment has been determined, then do it for real
	align($ref, $test, $offset[0], ($#{$test->{a_trace}} + $offset[0]));
	return 1;
}

sub align {
	# This takes the normalized traces and returns a best alignment of the two.
	# Rows are added to or deleted from the test trace to keep the alignment.
	# Best alignment is calculated by minimising the difference between the test
	# and reference sequence over the next 30 datapoints. It is adjusted every five
	# bases. Inserted lines are just a duplicate of the previous line.
	my ($ref, $test, $min_index, $trace_length) = @_;
	
	# Add/delete the appropriate number of lines to the test sequence to correct
	# for the offset value
	if ($min_index < 0) {
		# Need to add lines to test sequence, since it's behind
		for ($min_index..0) {
			for (@trace_list) {
				unshift @{$test->{$_}}, $test->{$_}[0];
			}
		}
	} elsif ($min_index > 0) {
		# Need to delete lines from test sequence, since it's ahead
		for (@trace_list) {
			splice @{$test->{$_}}, 0, $min_index;
		}
	}
	# Make a note of the offset value for datapoint numbering
	$test->{initial_offset} = $min_index;
	$ref->{initial_offset} = 0;
	
	
	# Now check alignments
	for my $index (0..($trace_length - 1)) {
		# Each third entry (starting from 1), check alignment
		next if $index%3;
		my $offset;
		my $start_pos = $index;
		my $end_pos = $index + 30;
		# Compare the scores in the current alignment with those one data point in 
		# either direction
		my $score = get_score($start_pos, $end_pos, 0, $ref, $test);
		my $pre_score = get_score($start_pos, $end_pos, -1, $ref, $test);
		my $post_score = get_score($start_pos, $end_pos, 1, $ref, $test);
		last if $score eq 'no_score' or $post_score eq 'no_score' or $pre_score eq 'no_score';
		# Work out offset
		# Default is 0; score is the lowest of the three
		if ($score < $pre_score && $score < $post_score) {$offset = 0}
		# if pre-score is the lowest, then offset = -1
		elsif ($pre_score < $score && $pre_score < $post_score) {$offset = -1}
		# if post-score is the lowest, then offset = 1
		elsif ($post_score < $score && $post_score < $pre_score) {$offset = 1}
		# If in doubt, default to no change
		else {$offset = 0}
		# Now insert or delete lines as required
		if ($offset == 1) {
		# The reference sample is behind, need to delete a row from test
			for (@trace_list) {
				splice @{$test->{$_}}, $index, 1;
			}
		}
		elsif ($offset == -1) {
		# The reference sample is ahead, need to add a row to test
			for (@trace_list) {
				splice @{$test->{$_}}, $index, 0, $test->{$_}[$index];
			}
		}
	}
	# Reset length of test sequence to correct value
	$test->{seq_length} = $#{$test->{a_trace}};
	return 1;
}
		
sub get_score {
	# Subroutine used in alignment testing - it gets the total difference between
	# the two submitted sections of array and returns it.
	my ($start, $end, $offset, $ref, $test) = @_;
	my $score;
	for my $index($start..$end) {
		return 'no_score' unless defined ($ref->{a_trace}[$index]) and defined($test->{a_trace}[$index + $offset]);
		for (@trace_list) {
			$score += abs($ref->{$_}[$index] - $test->{$_}[$index + $offset]);
		}
	}
	return $score;
}

sub differences {
	# Takes the two traces and calculates the difference between the two. Then 
	# squares this difference to highlight the larger (relevent) changes.
	# Also returns the length of the shortest sequence, so both are aligned on 
	# image generation
	my ($ref, $test) = @_;
	# Just compare up to end of shortest sequence
	my $min_index = $ref->{seq_length} < $test->{seq_length} ? $ref->{seq_length} : $test->{seq_length};
	 
	# Create a hash for the difference values
	my %diffs;
	# Loop through the traces and calculate the differences
	for my $index (0..$min_index) {
		# Need to do all four traces
		for (@trace_list) {
			# Get the difference
			my $diff = $ref->{$_}[$index] - $test->{$_}[$index];
			# Is the difference positive or negative (need to record, otherwise will 
			# be lost on squaring)
			my $sign = $diff < 0 ? -1 : 1;
			# Stop saturation by cutting off to a max value
			if (abs($diff) > 500) {$diff = 500}
			# Highlight differences by multiplying value by total values of OTHER
			# channels
			my $othersum = 1;
			my $samesum = 1;
			$diffs{$_}[$index] = $diff;
		}
		# Have now got difference for all four traces. Can accentuate real diffs
		for (@trace_list) {
			my $diff = $diffs{$_}[$index];
			my $sign = $diff < 0 ? -1 : 1;
			my $otherchannels = 1;
			# Sum all values in the other channels which have the opposite sign
			for my $channel (@trace_list) {
				next if $channel eq $_;
				my $value = $diffs{$channel}[$index];
				# Ignore if the sign is the same
				next if $value * $sign > 0;
				$otherchannels += $value;
			}
			my $finaldiff = ($sign * $diff * $diff * sqrt(abs($otherchannels))) / 5000;
			$finaldiff = $sign * 500 if abs($finaldiff) > 500;
			$diffs{$_}[$index] = $finaldiff;
		}
	}
	
	return $min_index, \%diffs;
}

sub get_image {
	# Routine to convert the supplied hash reference into a png image of a graph
	# This is the slowest section of the program by a long way (typically takes
	# around 75% of the total calculation time). The limiting factor appears to be
	# GD::Graph - optimization suggestions are welcome.
	my ($trace, $min, $max, $size, $output_file, $length, $chromatogram) = @_;
	
	my @labels;
	for (0..$length) {push @labels, ''}
	if ($chromatogram) {
		my @seq = split //, $trace->{sequence};
		for (0..@{$trace->{base_pos}}) {
			# Only mark every 10th base from 10.
			next if $_ % 10;
			next unless $_;
			last unless $trace->{base_pos}[$_];
			last if ($trace->{base_pos}[$_] - $trace->{initial_offset}) >= $length;
			# Correct for initial offset - puts numbers in same places as for 
			# chromatogram (otherwise they can get shifted)
			$labels[($trace->{base_pos}[$_-1] - $trace->{initial_offset})] = $_;
		}
	}
	
	my @plot_data = (\@labels, $trace->{a_trace}, $trace->{g_trace}, $trace->{c_trace}, $trace->{t_trace});
	
	# Define graph format
	my $mygraph = GD::Graph::lines->new($length, $size);
	$mygraph->set(
    x_label     => '',
		x_ticks     => 0,
	  y_label     => '',
		y_max_value => $max,
		y_min_value => $min,
		# Draw datasets in 'solid', 'dashed' and 'dotted-dashed' lines
		line_types  => [1, 1, 1, 1],
		# Set the thickness of line
		line_width  => 1,
    # Set colors for datasets
		dclrs       => ['green', 'black', 'blue', 'red'],
	) or warn $mygraph->error;
		
	my $myimage = $mygraph->plot(\@plot_data) or die $mygraph->error;
	open OUT, ">".$output_file;
	print OUT $myimage->png();
	close OUT;
	
	return 1;
}

sub demac {
	# Subroutine to strip out extra Mac headers
	for my $file ($temp_dir.$$."ref_seq", $temp_dir.$$."test_seq") {
		open IN, $file;
		open OUT, ">".$file."_new";
		while (<IN>) {
			if ($. == 1) {
				s/^.*ABIF/ABIF/;
				unless (/ABIF/) {error_end('abi_format')}
			}
			print OUT $_;
		}
		close IN;
		close OUT;
		rename $file."_new", $file
	}
}