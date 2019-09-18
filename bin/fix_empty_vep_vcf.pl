#!/usr/bin/perl -w
use strict;

my $header_file = "/data/bnf/ref/trusight_myeloid/vep_vcf_header";
my( $in, $out ) = ( $ARGV[0], $ARGV[1] );

unless( -s $in ) {
    print "Input VCF is empty...\n";
    exit;
}

unless( -s $out ) {

    my $vars = 0;
    open( IN, $in );
    while( <IN> ) {
	$vars++ unless /^#/;
    }
    close IN;

    if( $vars == 0 ) {
	print STDERR "Replacing empty output VCF with input VCF + VEP header.\n";
	if( $out =~ /\.gz$/ ) {
	    $out =~ s/\.gz$//;
	    system("head -n -1 $in > $out ; cat $header_file >> $out ; grep '^#CHROM' $in >> $out");
	    system("bgzip $out -f");
	}
	else {
	    system("head -n -1 $in > $out ; cat $header_file >> $out ; grep '^#CHROM' $in >> $out");
	}
    }
    else {
	print STDERR "WARNING: Input VCF has variants, but output VCF is empty. Something is fishy!\n";
    }
}


    
