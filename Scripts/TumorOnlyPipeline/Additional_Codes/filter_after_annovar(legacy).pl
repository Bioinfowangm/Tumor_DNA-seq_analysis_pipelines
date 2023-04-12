#!/usr/bin/perl 

use strict;
use warnings;

my $input = $ARGV[0];
my $output = $ARGV[1];
open I_f,$input;
open O_f,">$output";
Line: while(<I_f>){
    chomp;
    if (/Func.refGene/){
        print O_f "$_\n";
        next;
    }
    my @row = split /\t/,$_;

    # removing SNV/indels as Aaron suggested
    next if $row[5] eq 'intronic';
    next if $row[5] eq 'intergenic';
    next if $row[5] eq "ncRNA_exonic";
    next if $row[5] eq "ncRNA_intronic";
    next if $row[5] eq "UTR3";
    next if $row[5] eq "UTR5";
    next if $row[5] =~ /stream/ && $row[6] ne 'TERT';
    next if $row[8] eq 'synonymous SNV';

    # removing SNPs with >= 1% frequency in human populations
    for my $af(@row[10..28]){
        $af = 0 if $af eq '.';
        next unless $af;
        next Line if $af >=0.01;
    }

    # perform moderate filtering on ref/alt reads and MAF
    my ($dp,$ref,$alt);
    if($input =~ /FB/){
        my @info = split /:/,$row[-1];
        $dp = $info[2];
        $ref = $info[4];
        $alt = $info[6];
    }
    elsif($input =~ /HS/){
        my @dp4 = /DP4=(\d+),(\d+),(\d+),(\d+)/;
        $dp = $dp4[0]+$dp4[1]+$dp4[2]+$dp4[3];
        $ref = $dp4[0]+$dp4[1];
        $alt = $dp4[2]+$dp4[3];
    }
    elsif($input =~ /MT2/){
        my @info = split /:/,$row[-1];
        my @ad = split /,/,$info[1];
        $dp = $ad[0]+$ad[1];
        $ref = $ad[0];
        $alt = $ad[1];
    }
    next unless $dp >=20;
    next unless $alt/($alt+$ref) >= 0.03;

    print O_f "$_\n";
}
