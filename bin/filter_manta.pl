#!/usr/bin/perl -w


#
# Filter manta results with the follwoing filters
#  - Spanning paired-read support > 8% alt allele
#  - Split reads > 4% for the alt allele
#  - Min. depth 20
#


use strict;
use CMD::vcf;
use Data::Dumper;


my $vcf_to_parse = shift;

my( $meta, $var, $samples ) = CMD::vcf::parse_vcf($vcf_to_parse);


print $meta-> {'header_str'};

foreach my $sample (@{$samples}){
    
    foreach my $name (@{$var}){

	next if !($name ->{'GT'}->{$sample}->{'SR'});
	my @PR = split /,/, $name ->{'GT'}->{$sample}->{'PR'};
	my @SR = split /,/, $name ->{'GT'}->{$sample}->{'SR'};

	# minst 20 paired-läsningar, måste finnas både paired och split som stöttar SVn
	next if ($PR[1]+$PR[0] < 20 || $PR[1]+$PR[0]==0 || $SR[1]+$SR[0] == 0);

	# 
	if( $PR[1]/($PR[1]+$PR[0]) > .08 && $SR[1]/($SR[1]+$SR[0]) > .04) {
	    print $name->{'vcf_str'}."\n";
	}
	
    }
}



