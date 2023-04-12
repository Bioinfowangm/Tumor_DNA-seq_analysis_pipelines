#!/usr/bin/perl 

use strict;
use warnings;

my $input  = $ARGV[0];
my $output = $ARGV[1];
open I_f, $input;
open O_f, ">$output";

chomp(my @lines = <I_f>);

my $header = grep {/^\#CHROM/} @lines;
my @row  = split /\t/, $header;
if ( $row[-1] =~ /Normal/i ) {
    system("cp $input $output");
}
else {
    for my $l(@lines){
        my @row = split /\t/,$l;
        if (/^##/) {
            print O_f "$l\n";
        }
        else {
            print O_f
            join( "\t", @row[ 0 .. ( $#row - 2 ) ], $row[-1], $row[-2] ),
            "\n";
        }
    }
}
