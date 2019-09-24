#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use JSON;# qw( encode_json );

my $SID = $ARGV[0];
my $type = $ARGV[1];
my $coverage_file_summary = "cov_metrics.txt.sample_summary";
my $coverage_file = "cov_metrics.txt";
my $hs_metrics_file = "hs_metrics.txt";
my $align_metrics_file = "aln_metrics.txt";
my $insert_file = "is_metrics.txt";
my $dedup_metrics_file = $type."_metrics.txt";

my %results;
open( HS, $hs_metrics_file );
while( <HS> ) {
    if( /^\#SentieonCommandLine/ ) {
	    <HS>;
	    my $vals = <HS>;
	    my @a = split /\t/, $vals;
        $results{'pct_on_target'} = $a[18];
        #print "pct_on_target: $a[18]\n";
	    $results{'fold_enrichment'} = $a[25];
        #print "fold_enrichment: $a[25]\n";
    	$results{'mean_coverage'} = $a[22];
        #print "median_coverage: $a[22]\n";
    	$results{'fold_80'} = $a[32];
        #print "fold_80: $a[32]\n";
	}
}
close HS;

## INSSIZE ##
open( INS, $insert_file );
while( <INS> ) {
	if( /^\#SentieonCommandLine/ ) {
	    <INS>;
	    my $vals = <INS>;
	    my @a = split /\t/, $vals;
	    $results{'ins_size'} = $a[4];
        #print "ins_size: $a[4]\n";
	    $results{'ins_size_dev'} = $a[5];
        #print "ins_size_dev: $a[5]\n";
	}
}
close INS;


## DEDUP ##

open( DEDUP, $dedup_metrics_file );
while( <DEDUP> ) {
    if( /^\#SentieonCommandLine/ ) {
	    <DEDUP>;
	    my $vals = <DEDUP>;
	    my @a = split /\t/, $vals;
	    $results{'dup_reads'} = $a[6];
        #print "dup_reads: $a[6]\n";
	    $results{'num_reads'} = $a[2];
        #print "num_reads: $a[2]\n";
        $results{'dup_pct'} = $a[8];
        #print "dup_pct: $a[8]\n";
        my $mapped = $a[2]-$a[4];
        $results{'mapped_reads'} = $mapped;
        #print "mapped_reads: $mapped\n";
	}
}
close DEDUP;

my @cov;
open( COV, $coverage_file );
while( <COV> ) {

    unless( /^Locus/ ) {
        my $vals = <COV>;
        my @a = split /\t/, $vals;
        push @cov,$a[1];
    }
    
}
my @sorted_cov = sort @cov;

my $len_cov = scalar(@cov);
my $median;
if (($len_cov / 2) =~ /\d+\.\d+/) { 
    my $median_index = ($len_cov / 2) + 0.5;
    $median = $sorted_cov[$median_index-1];
}
else {
    my $median_index = ($len_cov / 2);
    $median = ($sorted_cov[$median_index-1] + $sorted_cov[$median_index-1]) / 2;
}
$results{'median_cov'} = $median;

## cov_% ##
my %pct_above_x;
open( COV_THRESH, $coverage_file_summary );
while( <COV_THRESH> ) {
    if( /^sample_id/ ) {
	    my $vals = <COV_THRESH>;
	    my @a = split /\t/, $vals;
        chomp @a;
        $pct_above_x{'1'} = $a[6];
        $pct_above_x{'10'} = $a[7];
        $pct_above_x{'30'} = $a[8];
        $pct_above_x{'100'} = $a[9];
        $pct_above_x{'250'} = $a[10];
        $pct_above_x{'500'} = $a[11];
	}
}
close COV_THRESH;

$results{'pct_above_x'} = \%pct_above_x;

my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );
