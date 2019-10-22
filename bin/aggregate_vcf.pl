#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CMD::vcf qw( parse_vcf );
use strict;
use Data::Dumper;



my %opt = ();
my %supported_callers = ( 'freebayes'=>1, 'mutect'=>1, 'vardict'=>1, 'pindel'=>1, 'gatkhc'=>1, 'manta'=>1, 'strelka'=>1, 'tnscope'=>1 );

GetOptions( \%opt, 'base=s', 'mutect=s@', 'freebayes=s@', 'vardict=s@', 'pindel=s@', 'manta=s@', 'gatkhc=s@', 'strelka=s@', 'tnscope=s@', 'excl-prefix=s', 'tumor-id=s', 'normal-id=s', 'fluffify-pindel' );
$opt{'excl-prefix'} = "DO_NOTHING" unless $opt{'excl-prefix'};

my( $base, $add ) = check_options( \%opt );


my @agg_vcf = aggregate_vcfs( $base, $add );

print join("\n", @agg_vcf);





sub aggregate_vcfs {
    my( $base, $add ) = @_;

    # Parse base VCF
    my( $base_meta, $base_vars, $base_samples ) = parse_vcf( $base->{file} );

    # Save base VCF's header in a string.
    my $base_header;

    if( is_gzipped($base->{file}) ) {
	open( BASE_VCF, "zcat ".$base->{file}."|" );
    }
    else {
	open( BASE_VCF, $base->{file} );
    }
    while( <BASE_VCF> ) {
	if( /^#/ ) {
	    if( /^##/ ) {
		$base_header .= $_;
	    }
	    else {

		# Change sample names in header for strelka
		if( $opt{base} eq "strelka" ) {
		    my @tmp = split /\t/;
		    $tmp[9] = $opt{'tumor-id'};
		    $tmp[10] = $opt{'normal-id'};
		    $_ = join( "\t", @tmp );

		    # Reverse order for TUMOR and NORMAL in strelka...
		    my $a = $base_samples->[0];
		    $base_samples->[0] = $base_samples->[1];
		    $base_samples->[1] = $a;
		}
		$base_header .= $_;
	    }
	}
    }
    chomp $base_header;

	
    

    # Make hash with real names. Needed because MuTect changes names to TUMOR and NORMAL and stores the real names in a meta field
    my( %base_names, %base_trans );
    foreach( @$base_samples ) {
	if( $base->{caller} eq "mutect" ) {
	    $base_names{ $_ } = &excl_prefix( $opt{'excl-prefix'}, $base_meta->{ SAMPLE }->{ $_ }->{ SampleName } );
	    $base_trans{ &excl_prefix( $opt{'excl-prefix'}, $base_meta->{ SAMPLE }->{ $_ }->{ SampleName } ) } = $_;
	}

	# Strelka requires that tumor and normal IDs are specified with the --tumor-id and --normal-id parameters!
	elsif( $base->{caller} eq "strelka" ) {
	    if( $_ eq "TUMOR" ) {
		$base_names{ $_ } = $opt{'tumor-id'};
		$base_trans{ $opt{'tumor-id'} } = $_;
	    }
	    elsif( $_ eq "NORMAL" ) {
		$base_names{ $_ } = $opt{'normal-id'};
		$base_trans{ $opt{'normal-id'} } = $_;
	    }
	    else {
		die "Strelka sample with ID other than NORMAL or TUMOR is not supported!"
	    }
	}	
	else {
	    $base_names{ $_ } = &excl_prefix( $opt{'excl-prefix'}, $_ );
	    $base_trans{ &excl_prefix( $opt{'excl-prefix'}, $_ ) } = $_;
	}
    }
 

  
    my @new;
    my %add_names;
    
    foreach my $add_data ( @$add ) {

	# Parse VCF to add to base VCF
	my( $add_meta, $add_vars, $add_samples ) = parse_vcf( $add_data->{file} );

	if( $opt{'fluffify-pindel'} and $add_data->{caller} eq "pindel" ) {
	    fluffify_pindel_variants( $add_vars );
	}
	
	# Make hash with real names. Needed because MuTect changes names to TUMOR and NORMAL. Stores the real names in a meta field
	foreach( @$add_samples ) {
	    if( $add_data->{caller} eq "mutect" and defined $add_meta->{ SAMPLE }->{ $_ }->{ SampleName } ) {
		$add_names{ $add_data->{file} }->{ &excl_prefix( $opt{'excl-prefix'}, $add_meta->{ SAMPLE }->{ $_ }->{ SampleName } ) } = $_;
	    }
	    elsif( $add_data->{caller} eq "strelka" ) {
		$add_names{ $add_data->{file} }->{ $opt{'normal-id'} } = $_ if $_ eq "NORMAL" and $opt{'normal-id'};
		$add_names{ $add_data->{file} }->{ $opt{'tumor-id'} } = $_ if $_ eq "TUMOR";
	    }
	    else {
		$add_names{ $add_data->{file} }->{ &excl_prefix( $opt{'excl-prefix'}, $_ ) } = $_;
	    }
	}

	for my $a_i ( 0..scalar(@$add_vars)-1 ) {
	    my $found = 0;

	    # Check if variant from 'add' is already found in 'base'
	    for my $b_i ( 0..scalar(@$base_vars)-1 ) {
		if( same_variant( $add_vars->[$a_i], $base_vars->[$b_i] ) ) {

		    # If the variant caller of the additional file is same as base, use add's GT data if better depth than base...
		    # If the additonal file has at least 3 times higher DP use this instead. FIXME: Handle this better? This is added because of Nextera XT in AML.
		    if( $add_vars->[$a_i]->{GT}->{$add_samples->[0]}->{DP} and $base_vars->[$b_i]->{GT}->{$base_samples->[0]}->{DP} and
			( $add_vars->[$a_i]->{GT}->{$add_samples->[0]}->{DP}  > 3 * $base_vars->[$b_i]->{GT}->{$base_samples->[0]}->{DP} ) ) {
			my %new_gt = %{ $add_vars->[$a_i]->{GT} };
			foreach my $sid ( %{ $add_vars->[$a_i]->{GT} } ) {
			    $new_gt{ excl_prefix( $opt{'excl-prefix'}, $sid ) } = $add_vars->[$a_i]->{GT}->{$sid};
			}
			$base_vars->[$b_i]->{GT} =  \%new_gt;
		    }

		    # Add 'add' caller info to FOUND_IN
		    push( @{ $base_vars->[$b_i]->{INFO}->{FOUND_IN} }, $add_data->{caller}."|".$add_data->{file} );
		    $found = 1;
		}
	    }

	    # If variant was not found in base set
	    unless( $found ) {

		# Check if it was previously found in an 'add' vcf
		my $new_found = 0;
		for my $n_i ( 0..scalar(@new)-1 ) {
		    if( same_variant( $add_vars->[$a_i], $new[$n_i] ) ) {
			
			# Add 'add' caller info to FOUND_IN
			push( @{ $new[ $n_i ]->{INFO}->{FOUND_IN} }, $add_data->{caller}."|".$add_data->{file} );
			$new_found = 1;
		    }
		}
		unless( $new_found ) {
		    push( @{ $add_vars->[$a_i]->{INFO}->{FOUND_IN} }, $add_data->{caller}."|".$add_data->{file} );
		    $add_vars->[$a_i]->{INFO}->{ORIGIN_FILE} = $add_data->{file};
		    $add_vars->[$a_i]->{INFO}->{ORIGIN_CALLER} = $add_data->{caller};
		    push( @new, $add_vars->[$a_i] );
		}
	    }
	}
    }

    
    my @final_vcf;
    push( @final_vcf, $base_header );
    # Output vcf lines from base vcf
    for my $b_i ( 0..scalar(@$base_vars)-1 ) {
	push( @{ $base_vars->[$b_i]->{INFO}->{FOUND_IN} }, $base->{caller}."|".$base->{file} );
	my $with_foundin = add_info_field( $base_vars->[$b_i]->{vcf_str}, "FOUND_IN=". join(",", @{ $base_vars->[$b_i]->{INFO}->{FOUND_IN} } ) );
#	print STDERR "WITH_FOUNDIN:".Dumper($with_foundin)."\nBASE SAMPLES:".Dumper($base_samples)."\nBASENAMES:".Dumper(\%base_names)."\nBASE_VARS".Dumper($base_vars->[$b_i])."\nBASETRANS".Dumper(\%base_trans)."\nBASE_CALLER". $base->{caller}."\n\n\n\n";
	my $final = reformat_GT( $with_foundin, $base_samples, \%base_names, $base_vars->[$b_i], \%base_trans, $base->{caller} );
	push( @final_vcf, $final );
    }

    # Output newly added vcf files from 'add' files
    for my $n_i ( 0..scalar(@new)-1 ) {
    	my $with_foundin = add_info_field( $new[$n_i]->{vcf_str}, "FOUND_IN=". join(",", @{ $new[$n_i]->{INFO}->{FOUND_IN} } ) );
	my $origin_file = $new[$n_i]->{INFO}->{ORIGIN_FILE};
	my $origin_caller = $new[$n_i]->{INFO}->{ORIGIN_CALLER};
	
	if( $origin_file =~ /germline/ ) {
	    $with_foundin = add_info_field( $with_foundin, "MYELOID_GERMLINE=1" );
	}

	my $final = reformat_GT( $with_foundin, $base_samples, \%base_names, $new[$n_i], $add_names{ $origin_file },  $origin_caller );

	push( @final_vcf, $final );
    }

    return @final_vcf;
}




sub reformat_GT {
    my( $vcf_str, $sample_order, $base_translate_names, $full_info, $translate_names, $caller ) = @_;

    my @vcf = split /\t/, $vcf_str;

    #print join( "\t", @vcf )."\n";
    my $out_vcf_str = join( "\t", @vcf[0..7] );
    $out_vcf_str .= "\tGT:DP:VD:AF";
    foreach my $sid ( @$sample_order ) {
	my $sample_name = $base_translate_names->{$sid};

	die "$sample_name not found in $caller vcf!". Dumper($translate_names) unless $translate_names->{ $sample_name };

	# MUTECT
	if( $caller eq "mutect" ) {

#	    print STDERR "--------------------\n".$caller."\n";
	    print STDERR Dumper($full_info) unless $full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
#	    print STDERR Dumper($translate_names);
#	    print STDERR Dumper($sample_name);

	    
	    
	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
	    $out_vcf_str .= ":".($ref_dp+$alt_dp);
	    $out_vcf_str .= ":".($alt_dp);
	    $out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{AF};
	}

	# FREEBAYES
	elsif( $caller eq "freebayes" ) {
	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
	    $out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
	    $out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{AO};
	    if( !$full_info->{GT}->{ $translate_names->{$sample_name} }->{AO} or $full_info->{GT}->{ $translate_names->{$sample_name} }->{AO} eq "." ) {
		$out_vcf_str .= ":0";
	    }
	    else {
		$out_vcf_str .= sprintf ":%.3f", $full_info->{GT}->{ $translate_names->{$sample_name} }->{AO} / $full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
	    }
	}

	# PINDEL
	elsif( $caller eq "pindel" ) {
	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
	    
	    $out_vcf_str .= ":".($ref_dp+$alt_dp);
	    $out_vcf_str .= ":".$alt_dp;
	    if( $alt_dp == 0 ) {
		$out_vcf_str .= ":0";
	    }
	    else {
		$out_vcf_str .= sprintf ":%.3f", $alt_dp / ($alt_dp + $ref_dp);
	    }
	}

	# MANTA
	elsif( $caller eq "manta" ) {
	    $out_vcf_str .= "\t0/1";
	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{PR};
	    $out_vcf_str .= ":".($ref_dp+$alt_dp);
            $out_vcf_str .= ":".$alt_dp;
            if( $alt_dp == 0 ) {
                $out_vcf_str .= ":0";
            }
            else {
                $out_vcf_str .= sprintf ":%.3f", $alt_dp / ($alt_dp + $ref_dp);
	    }	    
	}
	

	# GATK HaplotypeCaller
	elsif( $caller eq "gatkhc" ) {
	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
	    $out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
	    $out_vcf_str .= ":".$alt_dp;
	    $out_vcf_str .= ":".($alt_dp/$full_info->{GT}->{ $translate_names->{$sample_name} }->{DP});
	}
	    

	# STRELKA
	elsif( $caller eq "strelka" ) {

	    # INDELs
	    if( $full_info->{GT}->{ $translate_names->{$sample_name} }->{TAR} ) {
		my $ref_count = ( split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{TAR} )[0];
		my $alt_count = ( split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{TIR} )[0];
		my $dp = $full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
		
		my $af = 0;
		if( ($alt_count + $ref_count) > 0 ) {
		    $af = $alt_count / ( $alt_count + $ref_count );
		}
		my $gt = "0/1";
		$gt = "0/0" if $af < 0.01;
		$out_vcf_str .= "\t$gt:".$dp.":".$alt_count.":".$af;    
	    }

	    # SNVs
	    else {
#		print STDERR Dumper($full_info);
		my $REF_FIELD = $full_info->{REF}."U";
		my $ALT_FIELD = $full_info->{ALT}."U";
		#print STDERR "************* ".$full_info->{REF}."\t".$full_info->{ALT}."\t".$translate_names->{$sample_name}."\n";
		my $ref_count  = (split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{$REF_FIELD} )[0];
		my $alt_count  = (split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{$ALT_FIELD} )[0];
		my $dp = $full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
		my $af = $alt_count / ( $alt_count + $ref_count );
		my $gt = "0/1";
		$gt = "0/0" if $af < 0.01;		
		$out_vcf_str .= "\t$gt:".$dp.":".$alt_count.":".$af;    
	    }
	}
	
	# TNSCOPE
	elsif ( $caller eq "tnscope" ) {
		$out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
		my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
		$out_vcf_str .= ":".($ref_dp+$alt_dp);
		$out_vcf_str .= ":".($alt_dp);
		$out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{AF};

	}
	else {
	    die "GT-reformatting for $caller not yet implemented";
	}
    }

    return $out_vcf_str;
}

sub same_variant{
    my( $v1, $v2 ) = @_;

    print STDERR Dumper($v1)."\n.".Dumper($v2)."\n" unless $v1->{CHROM} and $v2->{CHROM} and $v1->{POS} and $v2->{POS} and $v1->{REF} and $v2->{REF} and $v1->{ALT} and $v2->{ALT};
    if( $v1->{CHROM} eq $v2->{CHROM} and
	$v1->{POS}   eq $v2->{POS} and
	$v1->{REF}   eq $v2->{REF} and
	$v1->{ALT}   eq $v2->{ALT} ) {

	return 1;
    }
    return 0;
}

sub check_options {
    my %opt = %{ $_[0] };

    help_text() unless $opt{base};
    
    my @files;
    my $base;
    my %seen_files;
    foreach my $opt_key ( sort keys %opt ) {
	if( $opt_key eq "base" ) {
	    unless( $supported_callers{ $opt{'base'} } ) {
		die "Base VCF set to '$opt{base}', which is not among the supported callers (". join(", ", keys %supported_callers).")"
	    }
	    unless( $opt{ $opt{base} } ) {
		die "Base VCF set to '$opt{base}', but no VCF was specified for that caller";
	    }
	}
	
	elsif ( $supported_callers{$opt_key} ) {
#	    print STDERR $_ ." ". join(",", keys %supported_callers)."\n";
	    my $first = 1;
	    foreach my $fn ( @{ $opt{$opt_key} } ) {
		die "$fn specified twice!" if $seen_files{ $fn };
		$seen_files{ $fn } = 1;

		unless( -s $fn ) {
		    print STDERR "WARNING: Variant file for '$opt_key', $fn could not be found.\n";# unless -s $fn;
		    next;
		}
		if( $opt{base} eq $opt_key and $first and $fn !~ /germline_myeloid/ ) {
		    $base = { 'file'=>$fn, 'caller'=>$opt_key };
		    $first = 0;
		}
		else {
		    push( @files, { 'file'=>$fn, 'caller'=>$opt_key } );
		}
		
	    }
	}
    }
    return ($base, \@files);
}


sub help_text {
    my $error = shift;
    
    print "\n\$ aggregate_vcf.pl --base BASE_VARIANT_CALLER [--freebayes/mutect/pindel/vardict/gatkhc/manta/strelka VCF [--fluffify-pindel] [--tumor-id --normal-id]\n\n";
    print "   --base        Defines which variant caller is used as base vcf (required)\n\n";

    print "Any number of these are allowed, at least one required:\n";
    print "   --freebayes VCF   Freebayes VCF path\n";
    print "   --vardict VCF     VarDict VCF path\n";
    print "   --mutect VCF      MuTect2 VCF path\n";
    print "   --manta VCF       Manta VCF path\n";
    print "   --pindel VCF      pindel VCF path\n";
    print "   --strelka VCF     Strelka2 VCF path\n";
	print "   --tnscope VCF     Tnscope VCF path\n";
    print "   --gatkhc VCF      GATK HaplotypeCaller VCF path\n\n";

    print "Additional options:\n";
    print "   --fluffify-pindel Modify pindel REF/ALT fields to be less than 1000 bp (comply with Manta).\n";
    print "   --tumor-id STR    REQUIRED IFF STRELKA: Sample ID for for tumor sample\n";
    print "   --normal-id STR   REQUIRED IFF STRELKA: Sample ID for for normal sample\n";
    print "\n";
    exit(0);
}

sub add_info_field {
    my $vcf_str = shift;
    my $data = shift;
    
    my @a = split /\t/, $vcf_str;
    $a[7] .= ";".$data;
    
    return join "\t", @a;
}


sub excl_prefix {
    my( $prefix, $name ) = @_;

    return $name if $prefix eq "DO_NOTHING";

    $name =~ s/^$prefix//;
    return $name;
}

sub is_gzipped {
    my $fn = shift;
    
    my $file_str = `file $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}


sub fluffify_pindel_variants {
    my $vars = shift;

    my $MAX_SIZE = 1000;
    
    foreach my $v ( @$vars ) {
	if( length($v->{REF}) >= $MAX_SIZE or length($v->{ALT}) >= $MAX_SIZE ) {
	    $v->{REF} = substr($v->{REF}, 0, 1);
	    $v->{ALT} = "<".$v->{INFO}->{SVTYPE}.">";
	}
    }
}
