#!/usr/bin/perl -w
use strict;
use CMD::vcf qw( parse_vcf );
use Data::Dumper;
use List::Util qw(sum);

my( $meta, $vars, $sample ) = parse_vcf( $ARGV[0] );

open( VCF, $ARGV[0] );
while( <VCF> ) {
    print $_ if $_ =~ /^#/;
}

#print Dumper($vars);exit;

my $T_ID = $sample->[0];

foreach my $var ( @$vars ) {

    next if $var->{FILTER} =~ /multi_event_alt_allele_in_normal/;
    my %T = %{ $var->{GT}->{$T_ID} };

    $T{DP} = sum( split /,/, $T{AD} );
    
    if( $T{DP} >= 20 and $T{AF} >= 0.04 ) {
	print $var->{ vcf_str }."\n";
    }
}




sub add_info_field {
    my $vcf_str = shift;
    my $data = shift;

    my @a = split /\t/, $vcf_str;
    $a[7] .= ";".$data;

    return join "\t", @a;
}
