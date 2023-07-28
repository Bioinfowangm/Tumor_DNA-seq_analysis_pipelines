#!/usr/bin/perl 

use strict;
use warnings;

my $input    = $ARGV[0];
my $output   = $ARGV[1];
my $ctl_name = $ARGV[2];
open I_f, $input;
open O_f, ">$output";

chomp( my @lines = <I_f> );

my $header = grep { /^\#CHROM/ } @lines;
my @row    = split /\t/, $header;
if ( $row[-1] eq $ctl_name ) {
    system("cp $input $output");
}
else {
    for my $l (@lines) {
        my @row = split /\t/, $l;
        if ( $l =~ /^##/ ) {
            print O_f "$l\n";
        }
        else {
            print O_f
              join( "\t", @row[ 0 .. ( $#row - 2 ) ], $row[-1], $row[-2] ),
              "\n";
        }
    }
}
