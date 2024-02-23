#!/usr/bin/perl
use strict;
use warnings;

# Check for correct number of arguments
die "Usage: $0 my.blastout output_file.txt

This script loads a BLAST output file and reports which OTUs still 'need' reblasting.
" unless @ARGV == 2;

# Initialize variables
my (%O, %G, %I, %L, @U);

my ($input_file, $output_file) = @ARGV;

# Open the input file
open(my $fh, '<', $input_file) or die "Could not open file '$input_file' $!";

# Open the output file
open(my $out, '>', $output_file) or die "Could not open file '$output_file' $!";

while (my $line = <$fh>) {
    chomp $line;
    next if $. < 2; # Skip the first line

    # Processing lines that match your criteria
    if ($line =~ /^(# Q|Otu)/) {
        if ($line !~ /^#/) {
            my @A = split /\t/, $line;
            next if exists $G{$A[1]};
            $O{$A[0]}++;
            $G{$A[1]}++;
            $I{$A[2]}++;
            $L{$A[4]}++;
        } elsif ($line =~ /^# Query/) {
            my @o = keys %O;
            my @i = keys %I;
            my @l = keys %L;
            if (@i == 1 || @l == 1) {
                print $out "$o[0]\n";
            }
            %O = ();
            %G = ();
            %I = ();
            %L = ();
        }
    }
}

close $fh;
close $out;
