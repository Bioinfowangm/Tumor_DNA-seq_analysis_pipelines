#!/usr/bin/perl 
#===========================================================================
#
#         FILE: checkorder.pl
#        USAGE: perl checkorder.pl
#
#       AUTHOR: Wang Meng, mengwang55@gmail.com
#      VERSION: 1.0
#      CREATED: 07/23/2019 09:31:31
#===========================================================================

use strict;
use warnings;

my $input  = $ARGV[0];
my $output = $ARGV[1];
open I_f, $input;
open O_f, ">$output";
my $line = `head -154 $input|tail -1`;
my @row  = split /\t/, $line;
if ( $row[-1] =~ /ormal/ ) {
    system("cp $input $output");
}
else {
    while (<I_f>) {
        chomp;
        my @row = split;
        if (/^##/) {
            print O_f "$_\n";
        }
        else {
            print O_f
              join( "\t", @row[ 0 .. ( $#row - 2 ) ], $row[-1], $row[-2] ),
              "\n";
        }
    }
}
