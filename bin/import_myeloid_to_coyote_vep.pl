#!/usr/bin/perl -w
use strict;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use CMD::vcf_arr qw( parse_vcf );
use Getopt::Long;
use JSON;

my %opt;
GetOptions( \%opt, 'vcf=s', 'id=s', 'clarity-sample-id=s', 'clarity-pool-id=s', 'bam=s', 'group=s', 'qc=s' );

my( $vcf, $id ) = ( $opt{vcf}, $opt{id} );
my @groups = split /,/, $opt{group};
#my @groups = ("mody");

# Read QC data
my @QC;
if( $opt{qc} ) {
    my @qc_files = split /,/, $opt{qc};
    foreach( @qc_files ) {
	if( -s $_ ) {
	    push @QC, read_json($_)
	}
	else {
	    print STDERR "WARNING: QC-json does not exist: $_\n.";
	}
    }
}


die "$vcf not found!" unless -s $vcf;

#################
# INSERT SAMPLE #
#################

# Connect to mongodb
my $client = MongoDB->connect();


# Prepare data to insert into sample collection
my $samples = $client->ns("coyote.samples");
my %sample_data = ( 'name'=>$id, 'groups'=>\@groups, 'time_added'=>DateTime->now, 'vcf_files'=>[$vcf] );
if ( scalar @QC > 0 ) {
    $sample_data{QC} = \@QC;
}

# Add clarity information if specified
if( $opt{'clarity-sample-id'} ) {
    $sample_data{'clarity-sample-id'} =  $opt{'clarity-sample-id'};
}
if( $opt{'clarity-pool-id'} ) {
    $sample_data{'clarity-pool-id'} =  $opt{'clarity-pool-id'};
}
if( $opt{'bam'} ) {
    $sample_data{'bam'} =  $opt{'bam'};
}

# Insert into collection
my $result2 = $samples->insert_one(\%sample_data);
my $SAMPLE_ID = $result2->inserted_id->value;
#print Dumper $result2;


print STDERR "ID:".$SAMPLE_ID."\n";


###################
# Insert variants #
###################
my( $meta, $data, $sample_order  ) = parse_vcf( $vcf ); #, $id );

# Fix/add various things to data structure
foreach( 0..scalar(@$data)-1 ) {

    # Add sample ID
    $data->[$_]->{SAMPLE_ID} = $SAMPLE_ID; #bless(\$SAMPLE_ID, "MongoDB::BSON::String");

    # Add type for pindel
    $data->[$_]->{INFO}->{TYPE} = $data->[$_]->{INFO}->{SVTYPE} if $data->[$_]->{INFO}->{SVTYPE};

    # Parse FOUND_IN
    my @arr;
    my @found_in = split /,/, $data->[$_]->{INFO}->{FOUND_IN};
    foreach( @found_in ) {
	my( $caller, $file ) = split /\|/;
	if( $file =~ /\/myeloid\// and $file =~ /germline/ ){
	    $caller .= "-germline";
	}
	elsif( $file =~ /\/myeloid\/N/ ){
	    $caller .= "-nextera";
	}
	elsif( $file =~/\/myeloid\/\d/ ) {
	    $caller .= "-panel";
	}
	push @arr, [ $caller, $file ];
    }
   
    $data->[$_]->{INFO}->{FOUND_IN} = \@arr;

    my $first = 1;
#    for my $sid ( @$sample_order ) {
    for my $i ( 0.. scalar( @{ $data->[$_]->{GT} } )-1 ) {
	
	my $gt = $data->[$_]->{GT}->[$i];

        if( !defined($gt->{AF}) or !defined($gt->{DP}) or !defined($gt->{VD}) or !defined($gt->{GT}) ) {
	    print STDERR Dumper($gt);
            die "Invalid VCF, should be aggregated with AF, DP, VD and GT";
        }

        if( $gt->{sample} =~ /^NORMAL_N/ ) {
            $gt->{sample} =~ s/^NORMAL_N//;
            $gt->{type} = "control";
        }
        elsif( $gt->{sample} =~ /^TUMOR_N/ ) {
            $gt->{sample} =~ s/^TUMOR_N//;
            $gt->{type} = "case";
        }
        else{
            $gt->{sample} =~ s/^N//;
            $gt->{type} = ( $first ? "case" : "control" );
        }
	$first = 0;
    }

    delete $data->[$_]->{vcf_str};
    delete $data->[$_]->{INFO}->{'technology.illumina'};
    
}    

#print Dumper($data);exit;


	 
 
my $variants = $client->ns("coyote.variants_idref");
my $var = $variants->with_codec( prefer_numeric => 1 );
my $result = $var->insert_many($data);
#print Dumper($result);


sub fix {
    my $str = shift;
    #$str =~ s/-/_/g;
    return $str;
}




sub read_json {
    my $fn = shift;

    print STDERR "Reading json $fn\n";

    open( JSON, $fn );
    my @json = <JSON>;
    my $decoded = decode_json( join("", @json ) );
    close JSON;

    return $decoded;
}
