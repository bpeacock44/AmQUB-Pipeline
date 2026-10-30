#!/usr/bin/env perl
use strict;
use warnings;
# Format 7, expected fields: qseqid sseqid pident length mismatch evalue bitscore ...
my ($input, $output, $limit) = @ARGV;
die "Usage: $0 input.blastout output.txt [max_target_seqs]\n" unless defined $input && defined $output;
$limit //= 1000;
open my $in, '<', $input or die "$input: $!";
open my $out, '>', $output or die "$output: $!";
my (%n, %best, %worst, @order);
while (my $line = <$in>) {
    next if $line =~ /^#/ || $line =~ /^\s*$/;
    chomp $line;
    my @f = split /\t/, $line;
    die "Invalid BLAST row (expected 10 fields): $line\n" unless @f == 10;
    my ($id, $score) = ($f[0], $f[6]);
    die "Invalid bitscore '$score' for $id\n" unless $score =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    push @order, $id unless exists $n{$id};
    ++$n{$id};
    $best{$id} = $score if !exists($best{$id}) || $score > $best{$id};
    $worst{$id} = $score if !exists($worst{$id}) || $score < $worst{$id};
}
for my $id (@order) {
    print {$out} "$id\n" if $n{$id} >= $limit && $worst{$id} == $best{$id};
}
close $out or die "Cannot close $output: $!";
