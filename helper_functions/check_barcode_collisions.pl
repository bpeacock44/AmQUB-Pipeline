#!/usr/bin/perl 

# Updated 12/5/23 by Beth and Paul.
# NOTE: Inline::C no longer needs to be installed to use this script.
# 

# NOTE: This script is step 1 of 3 scripts
# 1) check_barcode_collisions.pl
# 2) filter_barcode_noncollisions.py
# 3) fastq_convert_mm2pm_barcodes.py

use strict;
use Getopt::Std;
use lib '/helper_functions';

my @argv = @ARGV;
$0 =~ s/.*\/// if $0 =~ /\//;

#get options from command line
my $opt = {};
getopts('i:m:o:t:M:Cfh', $opt);

sub usage
{
	my $message = shift;
	print "
Determine usable and unusable mismatched barcodes in a fastq file (compared to barcodes in a mapping file)

Usage: $0 [options]
   
   [REQUIRED]
   -i <input file>   # The input fastq filepath (barcodes OR reads with barcodes in headers. See NOTE below)
   -m <map file>     # Barcodes in mapping file will be denoted in the output

   [OPTIONAL]
   -o <output file>  # Output filepath (default=uFQBC_<input_fastq_filename>_BC<barcode_length>_M<mismatches>.txt
   -M <integer>      # Mismatches to check (default=2)
   -t <integer>      # Terminate after this many barcodes have been loaded (for quick checks)
   -C                # Print mismatch Collisions (fastq barcodes with >= 1 bp mismatch to multiple mapping file barcodes)
   -f                # Force overwrite of output file
   -h                # This help

NOTE: Read fastq ID lines must end with a ':' and a barcode (eg., \":CTCGACTACTGA\")
IMPORTANT: The map file should contain the barcodes for ALL SAMPLES present in the fastq file or the results may contain undetected errors.

";
	exit;
}
usage() unless defined $$opt{i} && defined $$opt{m};
usage() if $$opt{h};
$$opt{M} = 2 unless $$opt{M};


$,=" ";

#get dirs of input files
my ($mapdir) = $$opt{m} =~ /(.*\/)/;
$mapdir = "./" unless $mapdir;
my ($fasqdir, $fastqfile) = $$opt{i} =~ /(.*\/)(.+)/;
$fasqdir = "./" unless($fasqdir);
$fastqfile = $$opt{i} unless $fastqfile;


#get barcodes from mapping file
my $map = get_mapping_file($$opt{m});
my $BCs = $$map{H}{"BarcodeSequence"};
my $bclen = length $$BCs[0];
my $mapBCcount = scalar @$BCs;


#create an output filename, saving in the same directory as the fastq file
$$opt{o} = "${fasqdir}uFQBC_${fastqfile}_BC${mapBCcount}_M$$opt{M}.txt" unless $$opt{o};
#die "** Warning ** Output file [$$opt{o}] already exists!\n" if -e $$opt{o} && ! $$opt{f};


#make a barcode->SampleID lookup
my %bcIDlu;
@bcIDlu{@$BCs} = @{$$map{H}{"#SampleID"}};
print "Mapping file [$$opt{m}] contains $mapBCcount barcodes\n";
print "Map barcode lengths=$bclen\n";


#determine input fastq type
open my $CHECK, "<", $$opt{i} or die "[$!]\n";
my $fqtype;
while(<$CHECK>) {
	if ($. % 4 == 2) {
		chomp;
		print "Seq_length=".length($_),"\n";
		if (length($_) > $bclen) {
			$fqtype = "sequence"
		} else {
			$fqtype = "barcode"
		}
	}
	last if $. > 1;
}
close $CHECK;
print "Fastq type: $fqtype\n";


#get .fastq barcodes
print "Reading fastq file [$$opt{i}]\n";
open my $FASTQ, "<", $$opt{i} or die "[$!]\n";
my %Fq;
if ($fqtype eq "sequence") {
	while(<$FASTQ>) {
		if ($. % 4 == 1) {
			#@M02457:252:000000000-BV354:1:1102:15647:1855 1:N:0:TATTCTGCGAGC
			chomp;
			/.+:(\w+)$/;
			$Fq{$1}++;
		}
		last if $$opt{t} && $. >= 4 * $$opt{t};
	}
} elsif ($fqtype eq "barcode") {
	while(<$FASTQ>) {
		if ($. % 4 == 2) {
			#TATTCTGCGAGC
			chomp;
			$Fq{$_}++;
		}
		last if $$opt{t} && $. >= 4 * $$opt{t};
	}
} else {

}
close $FASTQ;


#get a list of all barcodes found
my @FqBCs = sort keys %Fq;


#determine which FqBCs have mismatches (1 .. $$opt{M}) to any mapping file barcode (BCs)
print "Finding 0-$$opt{M} mismatches to mapfile barcodes\n";
my %C = %{find_mismatches($BCs, \@FqBCs)};


#save usable fastq barcodes to file
open SINGS, ">", $$opt{o} or die "[$!]\n";
$,=" ";

print "Checking for usable and unusable barcodes and saving results: ";
print SINGS "# $mapBCcount mapfile barcodes\n";
print SINGS "# Mapbc\tCount\t[SampleID]/\n";
print SINGS "# FQbc\tCount\tmmPos\n";

foreach my $bc (@$BCs)
{
	#initialize any mapfile bcs not found in the .fastq file
	$Fq{$bc} = 0 unless defined $Fq{$bc};
	
	local $,=" ";
	my $bccount;
	$bccount = $Fq{$bc};
	my $bcmmcount={};##
	
	#print the current mapping file barcode and its count in the fastq file
	print SINGS "\nbc\t$bc\t$Fq{$bc}\t[$bcIDlu{$bc}]\n";
	
	
	#check each mismatch (1bpmm, 2bpmm, ...) up to whatever $$opt{M} is
	foreach my $d (1 .. $$opt{M})
	{
		my $pmmm = "pm".$d."mm";
		
		#get a list of fastq barcodes with $d mismatches to the current mapping file barcode, and sorted... how?
		my @fqbcs = sortSub(\%{$C{$pmmm}{BC}{$bc}{FQ}}, \%Fq);#
		
		my $de = $d - 1;
		my @dees = (1 .. $de);
		
		#print the associated fastq barcodes (@fqbcs)
		foreach my $fqbc (@fqbcs)
		{
			#skip if it's a mapfile barcode (because those are handled by $bc)
			next if exists $C{PMC}{$fqbc};
			
			my @mapbcs = keys %{$C{$pmmm}{FQ}{$fqbc}{BC}}; #mapfile BCs > 1bp from this $fqbc?
			
			#determine if it's been printed
			my $alreadyprinted = 0;
			foreach my $dee(@dees)
			{
				$alreadyprinted = 1 if exists $C{multibc}{$fqbc} && exists $C{multibc}{$fqbc}{$de};
			}
			
			#skip if it's already been printed
			next if $alreadyprinted;
			
			
			#if it has more than one mismatch to more than one mapfile barcode:
			if ($$opt{C} && @mapbcs > 1)
			{
				#denote that this $fqbc can't be used by beginning the line with an asterisk
				print SINGS "*$d\t$fqbc\t$Fq{$fqbc}\t@{$C{$pmmm}{BC}{$bc}{FQ}{$fqbc}}";
				#print all of the [sampleIDs:bc] it maps to
				print SINGS "\t[ ";
				print SINGS "$bcIDlu{$_}:$_ " for @mapbcs;
				print SINGS "]\n";
			}
			#if it has no mismatch collisions
			if (@mapbcs == 1)
			{
				print SINGS "m$d\t$fqbc\t$Fq{$fqbc}\t@{$C{$pmmm}{BC}{$bc}{FQ}{$fqbc}}\n";
				$$bcmmcount{$d} += $Fq{$fqbc};
			}
		}
	}
	
	#tabulate counts of (non-colliding) mismatches
	my @bcmmtots;
	push @bcmmtots, $$bcmmcount{$_} for sort {$a<=>$b} keys %$bcmmcount;
	
	#to keep sample names lined up in 'tot' rows, add 0's whenever
	# no non-colliding mismatches exist below $$opt{M} mismatches
	my $cols_needed = $$opt{M} - (scalar @bcmmtots);
	push @bcmmtots, (0) x $cols_needed if $cols_needed;
	my $bcmmtots = join "\t", @bcmmtots;
	print SINGS "tot\t$bc\t$bccount\t$bcmmtots\t$bcIDlu{$bc}\n";
}
close SINGS;
print "[".$$opt{o}."]\n";


exit;



sub find_mismatches
{
	my ($BCs, $FqBCs) = @_;
	
	my %C;
	foreach my $bc (@$BCs)
	{
		#initialize
		$C{PMC}{$bc} = 0;
		
		#check all FqBCs with these regexes
		foreach my $fqbc (@$FqBCs)
		{
			my @diffpos = find_diffs($bc, $fqbc);
			
			#check if a perfect match
			if (!@diffpos)
			{
				$C{PMC}{$fqbc}++;  # just perfect matches
			}
			#check if 1bp or more mismatch
			elsif (@diffpos >= 1 && @diffpos <= $$opt{M})
			{
				my $d = @diffpos;
				my $pmmm = "pm".$d."mm";#create a 'key' for $d mismatches
				$C{$pmmm}{FQ}{$fqbc}{BC}{$bc} = \@diffpos;#map $fqbc->$bc
				$C{$pmmm}{BC}{$bc}{FQ}{$fqbc} = \@diffpos;#map $bc->$fqbc
				$C{multibc}{$fqbc}{$d}++;#keeps track of ANY fqbc that maps to more than 1 barcode with an arbitrary number of mismatches
			}
		}
	}
	return \%C;
}


sub sortSub
{
#$unsorted = \%{$C{$pmmm}{BC}{$bc}{FQ}} = $href{$fqbc}->\@diffpos
	my ($unsorted, $Fq) = @_;
	
	#create an array of indexes for the barcodes to be sorted
	my $bcs = keys %$unsorted;#the number of barcodes with $pmmm mismatches to the current barcode
	my @idxs2sort = (0 .. --$bcs);
	
	#sort
	my @sorted = sort
	{
		#
		for my $idx ( @idxs2sort )
		{
			#sort $fqbc's ascending by mismatch position, then descending by $fqbc counts
			my $cmp = $$unsorted{$a}->[$idx] <=> $$unsorted{$b}->[$idx]  ||  $Fq{$b} <=> $Fq{$a};
			return $cmp if $cmp;#returns 1 or -1 to 'sort' command?
		}
		return 0;
	} keys %$unsorted;
	
	return @sorted;
}


sub get_mapping_file
{
# load a mapping file and return a hashref with column headers
# as keys with values as refs to arrays of column values
# A partially-filled $map may be passed in and new map values
# will be added to it
# RETURNS:
#   $map{HEADERLIST} = \@headers;
#   @{$$map{H}{$header}}, $vals[$i];# Access column data by header name
#   $$map{ID}{$smplID} = \@vals;    # Access row data by sampleID
#   $$map{ColIDX}{$header} = $i;    # Index of $header in @headers

    my $file = shift;
    my $map = shift;
    $map = {} unless $map;
    
    open my $MAP, "<", $file or die "Can't open file '$file'! [$!]\n";
    
    #hack for newline difficulties, where $/ seems to have a chr 13 (CR) by default instead of a chr 10 (LF) that Linux uses
    local $/="\n";
    
    
    my $headers = <$MAP>;
    chomp $headers;
    my @headers = split /\t/, $headers;
    $$map{HEADERLIST} = \@headers;
    
    #my $r = 0;##
    while(my $line = <$MAP>)
    {
        next unless $line =~ /\S/;#skip blank lines
        chomp $line;
        my @vals = split( /\t/, $line );
        
        #access column data by column header names
        my $i = 0;
        foreach my $header(@headers)
        {
            push @{$$map{H}{$header}}, $vals[$i];
            $$map{ColIDX}{$header} = $i;
            $i++;
        }
        
        #access row data by array index
        #[This method makes column indexes match row indexes]
        my $smplID = $vals[0];
        $$map{ID}{$smplID} = \@vals;
        #$$map{RowIDX}{$smplID} = $r;##
        #$r++;##
        
        #TODO: check if other programs rely on this older method:
        #[It made column indexes ofset by +1 to row indexes]
        # #access row data by sampleID
        # my $smplID = shift @vals;
        # $$map{ID}{$smplID} = \@vals;
    }
    
    close $MAP;
    
    return $map;
}

sub find_diffs {
    my ($a, $b) = @_;
    my @diffs;

    my $length = length($a);
    
    for (my $i = 0; $i < $length; $i++) {
        if (substr($a, $i, 1) ne substr($b, $i, 1)) {
            push @diffs, $i;
        }
    }

    return @diffs;
}

#use Inline C => << 'EOC';
#//returns a reference to an array that contains mismatch positions (indexes)
#void find_diffs(char *a, char *b)
#{
#    int i, n;
#    Inline_Stack_Vars;
#    Inline_Stack_Reset;
#    for (i=0, n=0; *a && *b; a++, b++, i++)
#        if (*a != *b)
#            mXPUSHi(i), n++;
#    Inline_Stack_Return(n);
#}
#EOC
