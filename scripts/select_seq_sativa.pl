#!/usr/bin/perl -w

=head1 NAME

select_seq_sativa.pl - prepares input to, and processes results from, Sativa

=head1 USAGE

  perl select_seq_sativa.pl --stage <STAGE> --domain <DOMAIN>
    
  NOTE that names of input and output files need to be changed inside the script file itself!

=head1 REQUIRED ARGUMENTS

  --domain  <DOMAIN>    Specify domain. 'archaea' or 'bacteria'

  --stage  <STAGE>      Specify if the script is to prepare input to Sativa ('pre') or process results from Sativa ('post')

=head1 OPIONAL ARGUMENTS

  --help                Prints this help message

[Press q to close this help message]

=cut

use Getopt::Long;
use List::MoreUtils qw(uniq);

#$domain = "archaea";
#$domain = "bacteria";
$domain = undef;
$stage = undef;

&GetOptions('stage=s' => \$stage, 'domain=s' => \$domain, 'h!' => \$help);

if (!$stage) {
    system ('perldoc', $0);
    exit;
}
if (!$domain) {
    system ('perldoc', $0);
    exit;
}
if ($help) {
    system ('perldoc', $0);
    exit;
}

if ($domain ne "archaea") {
    if ($domain ne "bacteria") {
        print"Error: --domain needs to be set to either 'archaea' or 'bacteria'\n";
        exit;
    }
}

if ($stage eq "pre") {
    &make_sativa_input;
} elsif ($stage eq "post") {
    &make_sativa_filtered;
} else {
    print"Error: --stage needs to be set to either 'pre' or 'post'\n";
    exit;
}


#########################
### make sativa input ###

sub make_sativa_input {
$min_length = 1000;
$priotisetype = 1; # 1 is the standard mode
$numpertaxon = 5;
$max_rawlength = 2000;
$max_n = 0; # 0 is the standard mode
$file_midfix = "minRL".$max_rawlength."minAL".$min_length."maxPT".$numpertaxon;
## for archaea ##
if ($domain eq "archaea") {
    $gtdb_info_file = "gtdb_r202_16S/ar122_metadata_r202.tsv";
    $fwd_seq_file = "gtdb_r202_16S/archaea_ssu_all.fna";
    $fwd_alignment_file = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.alnfna";
    $rev_alignment_file = "gtdb_r202_16S/archaea_ssu_all.rev.hmmer.rfmask.alnfna";
    $seq_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.".$file_midfix.".alnfna";
    $tax_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.".$file_midfix.".tax";
    $log_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.".$file_midfix.".log";
}
## for bacteria ##
if ($domain eq "bacteria") {
    $gtdb_info_file = "gtdb_r202_16S/bac120_metadata_r202.tsv";
    $fwd_seq_file = "gtdb_r202_16S/bacteria_ssu_all.fna";
    $fwd_alignment_file = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.alnfna";
    $rev_alignment_file = "gtdb_r202_16S/bacteria_ssu_all.rev.hmmer.rfmask.alnfna";
    $seq_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.".$file_midfix.".alnfna";
    $tax_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.".$file_midfix.".tax";
    $log_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.".$file_midfix.".log";
}
&get_genome_info;
&read_alignment_fwd;
&read_alignment_rev;
&read_seq_fwd;
#&print_lengths_II;
&print_selection;
}

#################################
### print sativa-filtered set ###

sub make_sativa_filtered {
$skip_species = 1; # 1 is the standard mode
$priotisetype = 1; # 1 is the standard mode
$numpertaxon = 1; # 5 is the standard mode # OBS!

## for archaea ##
if ($domain eq "archaea") {
    $gtdb_info_file = "gtdb_r202_16S/ar122_metadata_r202.tsv";
    $fwd_alignment_file = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.alnfna";
    $rev_alignment_file = "gtdb_r202_16S/archaea_ssu_all.rev.hmmer.rfmask.alnfna";
    $sativa_infile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.minRL2000minAL1000maxPT5.alnfna";
    $sativa_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.minRL2000minAL1000maxPT5.mis";
    $fwd_seq_file = "gtdb_r202_16S/archaea_ssu_all.fna";
    $rev_seq_file = "gtdb_r202_16S/archaea_ssu_all.rev.fna";
    $seq_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.sativafilt".$numpertaxon."pt.fna";
    $log_outfile = "gtdb_r202_16S/archaea_ssu_all.hmmer.rfmask.sativafilt".$numpertaxon."pt.log";
}
## for bacteria ##
if ($domain eq "bacteria") {
    $gtdb_info_file = "gtdb_r202_16S/bac120_metadata_r202.tsv";
    $fwd_alignment_file = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.alnfna";
    $rev_alignment_file = "gtdb_r202_16S/bacteria_ssu_all.rev.hmmer.rfmask.alnfna";
    $sativa_infile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.minRL2000minAL1000maxPT5.alnfna";
    $sativa_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.minRL2000minAL1000maxPT5.mis";
    $fwd_seq_file = "gtdb_r202_16S/bacteria_ssu_all.fna";
    $rev_seq_file = "gtdb_r202_16S/bacteria_ssu_all.rev.fna";
    $seq_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.sativafilt".$numpertaxon."pt.fna";
    $log_outfile = "gtdb_r202_16S/bacteria_ssu_all.hmmer.rfmask.sativafilt".$numpertaxon."pt.log";
}
&get_genome_info;
&read_alignment_fwd;
&read_alignment_rev;
&read_sativa_input;
&read_sativa_output;
#&evaluate_sativa;
&print_sativa_filtered;
}
#############

sub get_genome_info {
    #print"Reading taxonomy\n";
    $num_types = 0;
    open (INFILE, $gtdb_info_file) || die ("Can't open $gtdb_info_file");
    while (<INFILE>) {
        $_ =~ s/\R//g;
        @fields = split(/\t/);
        #@subfields = split(/;/, $fields[1]);
        $id = $fields[0];
        $taxon = $fields[16];
        $type_strain = $fields[15];
        $id_taxon{$id} = $taxon;
        $taxon_n_id{$taxon}{$id} = 1;
        $id_type{$id} = $type_strain;
        $num_types++ if ($type_strain eq "t");
    }
}

sub read_seq_fwd {
    #print"Reading forward raw sequences\n";
    $num_seq = 0;
    open (INFILE, $fwd_seq_file) || die ("Can't open $fwd_seq_file");
    $seq = "";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            $num_seq++;
            if ($seq ne "") {
                $id_n_geneid_rawlen{$id}{$geneid} = length($seq);
                $id_n_geneid_ncount{$id}{$geneid} = $seq =~ tr/N//;
                #print "$id_n_geneid_ncount{$id}{$geneid}\n";
                #print "$id\t$geneid\t$id_type{$id}\t$id_n_geneid_rawlen{$id}{$geneid}\t$id_n_geneid_len{$id}{$geneid}\t$id_taxon{$id}\n";
                $seq = "";
            }
            @fields = split(/~/, $row); # split id at tilde
            $id = $fields[0];
            substr($id, 0, 1) = "";
            @fields = split(/\s+/, $fields[1]);
            $geneid = $fields[0];
        } else {
            $seq = $seq.$row;
        }
    }
    $id_n_geneid_rawlen{$id}{$geneid} = length($seq);
    $id_n_geneid_ncount{$id}{$geneid} = $seq =~ tr/N//;
    $seq = "";
    close (INFILE);
}

sub read_alignment_fwd {
    #print"Reading forward aligned sequences\n";
    open (INFILE, $fwd_alignment_file) || die ("Can't open $fwd_alignment_file");
    $seq = "";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            if ($seq ne "") {
                $id_n_geneid_seq{$id}{$geneid} = $seq;
                $seq_mod = $seq;
                $seq_mod =~ tr/-//d;
                #print"\n$seq\n";
                $id_n_geneid_len{$id}{$geneid} = length($seq_mod);
                $id_n_geneid_dir{$id}{$geneid} = "fwd";
                $seq = "";
            }
            @fields = split(/~/, $row); # split id at tilde
            $id = $fields[0];
            substr($id, 0, 1) = "";
            @fields = split(/\s+/, $fields[1]);
            $geneid = $fields[0];
            #print"$id\n"; die;
        } else {
            $seq = $seq.$row;
        }
    }
    $id_n_geneid_seq{$id}{$geneid} = $seq;
    $seq_mod = $seq;
    $seq_mod =~ tr/-//d;
    $id_n_geneid_len{$id}{$geneid} = length($seq_mod);
    $id_n_geneid_dir{$id}{$geneid} = "fwd";
    $seq = "";
    close (INFILE);
}

sub read_alignment_rev {
    #print"Reading reverse aligned sequences\n";
    open (INFILE, $rev_alignment_file) || die ("Can't open $rev_alignment_file");
    $seq = "";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            if ($seq ne "") {
                $seq_mod = $seq;
                $seq_mod =~ tr/-//d;
                #print"\n$seq\n";
                $len = length($seq_mod);
                if ($len > $id_n_geneid_len{$id}{$geneid}) {
                    $id_n_geneid_seq{$id}{$geneid} = $seq;
                    $id_n_geneid_len{$id}{$geneid} = $len;
                    $id_n_geneid_dir{$id}{$geneid} = "rev";
                }
                $seq = "";
            }
            @fields = split(/~/, $row); # split id at tilde
            $id = $fields[0];
            substr($id, 0, 1) = "";
            @fields = split(/\s+/, $fields[1]);
            $geneid = $fields[0];
            #print"$id\n"; die;
        } else {
            $seq = $seq.$row;
        }
    }
    $seq_mod = $seq;
    $seq_mod =~ tr/-//d;
    #print"\n$seq\n";
    $len = length($seq_mod);
    if ($len > $id_n_geneid_len{$id}{$geneid}) {
        $id_n_geneid_seq{$id}{$geneid} = $seq;
        $id_n_geneid_len{$id}{$geneid} = $len;
        $id_n_geneid_dir{$id}{$geneid} = "rev";
    }
    $seq = "";
    close (INFILE);
}

sub read_sativa_input {
    open (INFILE, $sativa_infile) || die ("Can't open $sativa_infile");
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            substr($row, 0, 1) = "";
            @fields = split(/\|/, $row);
            $id_geneid_included{$fields[0]} = $fields[1];
        }
    }
}

sub read_sativa_output {
    open (INFILE, $sativa_outfile) || die ("Can't open $sativa_outfile");
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        next if (substr($row, 0, 1) eq ";");
        @fields = split(/\s+/, $row);
        @subfields = split(/\|/, $fields[0]);
        if ($skip_species == 1) {
            next if ($fields[1] eq "Species");
        }
        $id_geneid_mislab{$subfields[0]} = $subfields[1];
    }
}

sub print_lengths_I {
    foreach $id (keys %id_n_geneid_seq) {
        $longest_geneid = "N/A";
        $longest_length = -1;
        foreach $geneid (keys %{$id_n_geneid_seq{$id}}) {
            if ($id_n_geneid_len{$id}{$geneid} > $longest_length) {
                $longest_length = $id_n_geneid_len{$id}{$geneid};
                $longest_geneid = $geneid;
                #print"$id\t$geneid\t$id_n_geneid_len{$id}{$geneid}\t$id_taxon{$id}\t$id_type{$id}\n";
            }
        }
        next if (!defined $id_type{$id}); # if seq gene ID absent in gtdb metadata file
        print"$id\t$longest_geneid\t$id_n_geneid_len{$id}{$longest_geneid}\t$id_type{$id}\t$id_taxon{$id}\n";
        #print"$id\t$longest_geneid\t$id_n_geneid_len{$id}{$longest_geneid}\t$id_type{$id}\n";
    }
}

sub print_lengths_II {
    foreach $id (keys %id_n_geneid_seq) {
        $longest_geneid = "N/A";
        $longest_length = -1;
        foreach $geneid (keys %{$id_n_geneid_seq{$id}}) {
            if ($id_n_geneid_len{$id}{$geneid} > $longest_length) {
                $longest_length = $id_n_geneid_len{$id}{$geneid};
                $longest_geneid = $geneid;
                #print"$id\t$geneid\t$id_n_geneid_len{$id}{$geneid}\t$id_taxon{$id}\t$id_type{$id}\n";
            }
        }
        next if (!defined $id_type{$id}); # if seq gene ID absent in gtdb metadata file
        $length_n_taxon{$longest_length}{$id_taxon{$id}} = 1;
    }
    $num_taxa = 0;
    for ($i = 1; $i < $min_length; $i++) {
        $num_taxa = 0;
        if (defined $length_n_taxon{$i}) {
            $num_taxa = (keys %{$length_n_taxon{$i}});
        }
        print"$i\t$num_taxa\n";
            #foreach $taxon (keys %{$length_n_taxon{$i}}) {
            #    $sofartaxon{$taxon} = 1;
            #}
            #$num_taxa = (keys %sofartaxon);
        #}
        #print"$i\t$num_taxa\n";
    }
}

sub print_selection {
    open (OUT1, ">$seq_outfile");
    open (OUT2, ">$tax_outfile");
    open (OUT3, ">$log_outfile");
    $num_taxa = 0;
    $num_printed_taxa = 0;
    $num_printed_seqs = 0;
    $num_printed_types = 0;
    foreach $taxon (keys %taxon_n_id) {
        $num_taxa++;
        #print"$taxon\n";
        @ids = (keys %{$taxon_n_id{$taxon}});
        %id_longestlen = ();
        %id_longestgeneid = ();
        $typeid = undef;
        foreach $id (@ids) {
            #die if (!defined $id_type{$id}); # should not happen
            $longest_geneid = "N/A";
            $longest_length = -1;
            #$fewest_n = 1000000;
            
            # rank by length
            foreach $geneid (keys %{$id_n_geneid_seq{$id}}) {
                if ($id_n_geneid_len{$id}{$geneid} > $longest_length) {
                    next if ($id_n_geneid_rawlen{$id}{$geneid} > $max_rawlength);
                    next if ($id_n_geneid_ncount{$id}{$geneid} > $max_n);
                    $longest_length = $id_n_geneid_len{$id}{$geneid};
                    $longest_geneid = $geneid;
                }
            }
            next if ($longest_length < $min_length);
            
            # rank by nu of N and length
            #foreach $geneid (keys %{$id_n_geneid_seq{$id}}) {
            #    next if ($id_n_geneid_rawlen{$id}{$geneid} > $max_rawlength);
            #    next if ($id_n_geneid_len{$id}{$geneid} < $min_length);
            #    if ($id_n_geneid_ncount{$id}{$geneid} < $fewest_n) {
            #        $fewest_n = $id_n_geneid_ncount{$id}{$geneid};
            #        $longest_length = $id_n_geneid_len{$id}{$geneid};
            #        $longest_geneid = $geneid;
            #    } elsif $id_n_geneid_ncount{$id}{$geneid} == $fewest_n) {
            #        if ($id_n_geneid_len{$id}{$geneid} > $longest_length) {
            #            $longest_length = $id_n_geneid_len{$id}{$geneid};
            #            $longest_geneid = $geneid;
            #        }
            #    }
            #}

            $id_longestlen{$id} = $longest_length;
            $id_longestgeneid{$id} = $longest_geneid;
            if ($id_type{$id} eq "t") {
                $typeid = $id;
                if ($priotisetype == 1) {
                    $id_longestlen{$id} = 1000000;
                    $num_printed_types++;
                }
            }
            
        }
        @sortids = (sort { $id_longestlen{$b} <=> $id_longestlen{$a} } keys %id_longestlen);
        $num_printed_taxa++ if (@sortids > 0);
        for ($i = 0; $i < @sortids; $i++) {
            next if ($i >= $numpertaxon);
            $id = $sortids[$i];
            print OUT1 ">$id|$id_longestgeneid{$id}\n$id_n_geneid_seq{$id}{$id_longestgeneid{$id}}\n";
            print OUT2 "$id|$id_longestgeneid{$id}\t$id_taxon{$id}\n";
            #print "$id|$id_longestgeneid{$id}\t$id_taxon{$id}\t$id_n_geneid_rawlen{$id}{$id_longestgeneid{$id}}\n";
            $num_printed_seqs++;
        }
    }
    print OUT3 "\nInput GTDB info file: $gtdb_info_file\n";
    print OUT3 "Input Raw sequence file: $fwd_seq_file\n";
    print OUT3 "Input Fwd-aligned sequence file: $fwd_alignment_file\n";
    print OUT3 "Input Rev-aligned sequence file: $rev_alignment_file\n";
    print OUT3 "Output aligned sequence file: $seq_outfile\n";
    print OUT3 "Output taxonomy file: $tax_outfile\n";
    print OUT3 "\nMax raw length was set to: $max_rawlength\n";
    print OUT3 "Max number of N was set to: $max_n\n";
    print OUT3 "Min aligned length was set to: $min_length\n";
    print OUT3 "Max included sequences per taxon was set to: $numpertaxon\n";
    print OUT3 "Prioritise type sequences was set to: $priotisetype\n";
    print OUT3 "\nTotal number of taxa: $num_taxa\n";
    print OUT3 "Total number of sequences: $num_seq\n";
    print OUT3 "Total number of type genomes: $num_types\n";
    print OUT3 "Filtered number of taxa: $num_printed_taxa\n";
    print OUT3 "Filtered number of sequences: $num_printed_seqs\n";
    print OUT3 "Filtered number of type sequences: $num_printed_types\n";
    
    print "\nInput GTDB info file: $gtdb_info_file\n";
    print "Input Raw sequence file: $fwd_seq_file\n";
    print "Input Fwd-aligned sequence file: $fwd_alignment_file\n";
    print "Input Rev-aligned sequence file: $rev_alignment_file\n";
    print "Output aligned sequence file: $seq_outfile\n";
    print "Output taxonomy file: $tax_outfile\n";
    print "\nMax raw length was set to: $max_rawlength\n";
    print "Max number of N was set to: $max_n\n";
    print "Min aligned length was set to: $min_length\n";
    print "Max included sequences per taxon was set to: $numpertaxon\n";
    print "Prioritise type sequences was set to: $priotisetype\n";
    print "\nTotal number of taxa: $num_taxa\n";
    print "Total number of sequences: $num_seq\n";
    print "Total number of type genomes: $num_types\n";
    print "Filtered number of taxa: $num_printed_taxa\n";
    print "Filtered number of sequences: $num_printed_seqs\n";
    print "Filtered number of type sequences: $num_printed_types\n";
    close(OUT1);
    close(OUT2);
    close(OUT3);
}

sub evaluate_sativa {
    foreach $id (keys %id_n_geneid_seq) {
        $all_taxa{$id_taxon{$id}} = 1;
    }
    foreach $id (keys %id_geneid_included) {
        $sativa_input_taxa{$id_taxon{$id}} = 1;
        if (!defined $id_geneid_mislab{$id}) {
            $sativa_output_taxa{$id_taxon{$id}} = 1;
        }
    }
    $num_all = (keys %all_taxa);
    $num_sativa_input = (keys %sativa_input_taxa);
    $num_sativa_output = (keys %sativa_output_taxa);
    print"\nNumber of taxa total: $num_all\nNumber of taxa pre-sativa: $num_sativa_input\nNumber of taxa post-sativa: $num_sativa_output\n";

    $num_sativa_input = (keys %id_geneid_included);
    $num_sativa_output = $num_sativa_input - (keys %id_geneid_mislab);
    print"\nNumber of sequences pre-sativa: $num_sativa_input\nNumber of sequences post-sativa: $num_sativa_output\n";
    
    foreach $taxon (keys %sativa_input_taxa) {
        if (!defined $sativa_output_taxa{$taxon}) {
            push(@missed_taxa, $taxon);
        }
    }
    @missed_taxa = sort(@missed_taxa);
    #print"\nTaxa lost in Sativa:\n";
    foreach $taxon (@missed_taxa) {
    #    print"$taxon\n";
    }
}

sub print_sativa_filtered {
    $num_taxa = 0;
    $num_printed_taxa = 0;
    $num_printed_seqs = 0;
    $num_printed_types = 0;
    foreach $taxon (keys %taxon_n_id) {
        $num_taxa++;
        @ids = (keys %{$taxon_n_id{$taxon}});
        %id_len = ();
        $typeid = undef;
        foreach $id (@ids) {
            next if (!defined $id_geneid_included{$id});
            next if (defined $id_geneid_mislab{$id});
            $geneid = $id_geneid_included{$id};
            $id_len{$id} = $id_n_geneid_len{$id}{$geneid};
            if ($id_type{$id} eq "t") {
                $typeid = $id;
                if ($priotisetype == 1) {
                    $id_len{$id} = 1000000;
                }
            }
        }
        @sortids = (sort { $id_len{$b} <=> $id_len{$a} } keys %id_len);
        $num_printed_taxa++ if (@sortids > 0);
        for ($i = 0; $i < @sortids; $i++) {
            next if ($i >= $numpertaxon);
            $id = $sortids[$i];
            $id_geneid_filtered{$id} = $id_geneid_included{$id};
            $num_printed_seqs++;
            $num_printed_types++ if ($id_type{$id} eq "t");
        }
    }
    
    open (OUT1, ">$seq_outfile");
    # printing filtered fwd seqs
    open (INFILE, $fwd_seq_file) || die ("Can't open $fwd_seq_file");
    $seq = "";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            if ($seq ne "") {
                if (defined $id_geneid_filtered{$id}) {
                    if ($id_geneid_filtered{$id} eq $geneid) {
                        if ($id_n_geneid_dir{$id}{$geneid} eq "fwd") {
                            print OUT1 ">$name\n$seq\n"
                        }
                    }
                }
                $seq = "";
            }
            substr($row, 0, 1) = "";
            $name = $row;
            @fields = split(/~/, $row); # split id at tilde
            $id = $fields[0];
            @fields = split(/\s+/, $fields[1]);
            $geneid = $fields[0];
        } else {
            $seq = $seq.$row;
        }
    }
    if (defined $id_geneid_filtered{$id}) {
        if ($id_geneid_filtered{$id} eq $geneid) {
            if ($id_n_geneid_dir{$id}{$geneid} eq "fwd") {
                print OUT1 ">$name\n$seq\n"
            }
        }
    }
    $seq = "";
    close (INFILE);
    
    # printing filtered rev seqs
    open (INFILE, $rev_seq_file) || die ("Can't open $rev_seq_file");
    $seq = "";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            if ($seq ne "") {
                if (defined $id_geneid_filtered{$id}) {
                    if ($id_geneid_filtered{$id} eq $geneid) {
                        if ($id_n_geneid_dir{$id}{$geneid} eq "rev") {
                            print OUT1 ">$name\n$seq\n"
                        }
                    }
                }
                $seq = "";
            }
            substr($row, 0, 1) = "";
            $name = $row;
            @fields = split(/~/, $row); # split id at tilde
            $id = $fields[0];
            @fields = split(/\s+/, $fields[1]);
            $geneid = $fields[0];
        } else {
            $seq = $seq.$row;
        }
    }
    if (defined $id_geneid_filtered{$id}) {
        if ($id_geneid_filtered{$id} eq $geneid) {
            if ($id_n_geneid_dir{$id}{$geneid} eq "rev") {
                print OUT1 ">$name\n$seq\n"
            }
        }
    }
    $seq = "";
    close (INFILE);
    close(OUT1);
    
    open (OUT2, ">$log_outfile");
    print OUT2 "\nInput GTDB info file: $gtdb_info_file\n";
    print OUT2 "Input Fwd raw sequence file: $fwd_seq_file\n";
    print OUT2 "Input Rev raw sequence file: $rev_seq_file\n";
    print OUT2 "Input sativa infile: $sativa_infile\n";
    print OUT2 "Input sativa outfile: $sativa_outfile\n";
    print OUT2 "Output sequence file: $seq_outfile\n";
    print OUT2 "\nMax sequences per taxon was set to: $numpertaxon\n";
    print OUT2 "Prioritise type sequences was set to: $priotisetype\n";
    print OUT2 "Ignore misclassifications at species level was set to: $skip_species\n";
    print OUT2 "\nTotal number of taxa: $num_taxa\n";
    print OUT2 "Filtered number of taxa: $num_printed_taxa\n";
    print OUT2 "Filtered number of sequences: $num_printed_seqs\n";
    print OUT2 "Filtered number of type sequences: $num_printed_types\n";
    close(OUT2);
    
    print "\nInput GTDB info file: $gtdb_info_file\n";
    print "Input Fwd raw sequence file: $fwd_seq_file\n";
    print "Input Rev raw sequence file: $rev_seq_file\n";
    print "Input sativa infile: $sativa_infile\n";
    print "Input sativa outfile: $sativa_outfile\n";
    print "Output sequence file: $seq_outfile\n";
    print "\nMax sequences per taxon was set to: $numpertaxon\n";
    print "Prioritise type sequences was set to: $priotisetype\n";
    print "Ignore misclassifications at species level was set to: $skip_species\n";
    print "\nTotal number of taxa: $num_taxa\n";
    print "Filtered number of taxa: $num_printed_taxa\n";
    print "Filtered number of sequences: $num_printed_seqs\n";
    print "Filtered number of type sequences: $num_printed_types\n";
}


######## end ########

