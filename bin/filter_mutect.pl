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
my $N_ID = $sample->[1];

foreach my $var ( @$vars ) {
    #print Dumper($var)."\n";
#    next if $var->{FILTER} =~ /multi_event_alt_allele_in_normal/;
    my %T = %{ $var->{GT}->{$T_ID} };
    my %N = %{ $var->{GT}->{$N_ID} };

    $T{DP} = sum( split /,/, $T{AD} );
    $N{DP} = sum( split /,/, $N{AD} );
    
    if( $T{DP} >= 500 and $N{DP} >= 100
	and $T{AF} >= 0.05 and $N{AF} <= 0.02
	and $var->{INFO}->{'Func.refGene'} ne "intronic" ) {

	print add_info_field( $var->{ vcf_str }, "MUTECT_STATUS=good" )."\n";
    }
    elsif( $T{DP} * $T{AF} >= 12 and $N{AF} * 3 < $T{AF} and $T{AF} > 0.025 and $N{DP} >= 10 ) {
	print add_info_field( $var->{ vcf_str }, "MUTECT_STATUS=maybe" )."\n";
    }
}




sub add_info_field {
    my $vcf_str = shift;
    my $data = shift;

    my @a = split /\t/, $vcf_str;
    $a[7] .= ";".$data;

    return join "\t", @a;
}
