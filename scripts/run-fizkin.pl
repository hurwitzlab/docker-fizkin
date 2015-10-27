#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use autodie;
use Algorithm::Numerical::Sample 'sample';
use Bio::SeqIO;
use Data::Dump 'dump';
use File::Basename qw'dirname basename';
use File::Copy;
use File::Find::Rule;
use File::Path 'make_path';
use File::RandomLine;
use File::Spec::Functions;
use File::Temp 'tempfile';
use Getopt::Long;
use List::MoreUtils 'uniq';
use Math::Combinatorics 'combine';
use Pod::Usage;
use Readonly;
use Statistics::Descriptive::Discrete;

my $DEBUG = 0;

main();

# --------------------------------------------------
sub main {
    my %args = get_args();

    if ($args{'help'} || $args{'man_page'}) {
        pod2usage({
            -exitval => 0,
            -verbose => $args{'man_page'} ? 2 : 1
        });
    }

    unless ($args{'in_dir'}) {
        pod2usage('No input directory');
    }

    unless ($args{'out_dir'}) {
        pod2usage('No output directory');
    }

    unless (-d $args{'in_dir'}) {
        pod2usage("Bad input dir ($args{'in_dir'})");
    }

    unless (-d $args{'out_dir'}) {
        make_path($args{'out_dir'});    
    }

    if ($args{'debug'}) {
        $DEBUG = 1;
    }

    run(\%args);
}

# --------------------------------------------------
sub debug {
    if ($DEBUG && @_) {
        say @_;
    }
}

# --------------------------------------------------
sub get_args {
    my %args = (
        kmer_size   => 20,
        mode_min    => 1,    
        num_threads => 12,
        hash_size   => '100M',
        max_seqs    => 300_000,
        max_samples => 15,
    );

    GetOptions(
        \%args,
        'in_dir=s',
        'out_dir=s',
        'kmer_size:i',
        'mode_min:i',
        'num_threads:i',
        'hash_size:s',
        'max_seqs:i',
        'max_samples:i',
        'help',
        'man',
    ) or pod2usage(2);

    return %args;
}

# --------------------------------------------------
sub sys_exec {
    my @args = @_;

    debug("exec = ", join(' ', @args));

    unless (system(@args) == 0) {
        die sprintf("Cannot execute %s: %s", join(' ', @args), $?);
    }

    return 1;
}

# --------------------------------------------------
sub run {
    my $args = shift;

    say "Subsetting";
    subset_files($args);

    say "Indexing";
    jellyfish_index($args);

    say "Kmerize";
    kmerize($args);

    say "Pairwise comp";
    pairwise_cmp($args);

    say "Matrix";
    make_matrix($args);
}

# --------------------------------------------------
sub jellyfish_index {
    my $args       = shift;
    my $subset_dir = catdir($args->{'out_dir'}, 'subset');
    my $jf_idx_dir = catdir($args->{'out_dir'}, 'jf');

    unless (-d $jf_idx_dir) {
        make_path($jf_idx_dir);
    }

    my $mer_size  = $args->{'kmer_size'};
    my $hash_size = $args->{'hash_size'};
    my $threads   = $args->{'num_threads'};

    my @files = File::Find::Rule->file()->in($subset_dir);

    printf "Found %s files in %s\n", scalar(@files), $subset_dir;

    my $file_num = 0;
    for my $file (@files) {
        my $basename = basename($file);
        printf "%5d: %s, ", ++$file_num, $basename;

        my $jf_file = catfile($jf_idx_dir, $basename . '.jf');

        if (-e $jf_file) {
            say "index exists.";
        }
        else {
            say "indexing";
            sys_exec('jellyfish', 'count', 
                '-m', $mer_size, 
                '-s', $hash_size,
                '-t', $threads,
                '-o', $jf_file,
                $file
            );
        }
    }
}

# --------------------------------------------------
sub kmerize {
    my $args       = shift;
    my $subset_dir = catdir($args->{'out_dir'}, 'subset');
    my $kmer_dir   = catdir($args->{'out_dir'}, 'kmer');

    unless (-d $kmer_dir) {
        make_path($kmer_dir);
    }

    my $mer_size = $args->{'kmer_size'};
    my @files    = File::Find::Rule->file()->in($subset_dir);

    printf "Found %s files in %s\n", scalar(@files), $subset_dir;

    my $file_num = 0;
    FILE:
    for my $file (@files) {
        my $basename = basename($file);
        printf "%5d: %s, ", ++$file_num, $basename;

        my $kmer_file = catfile($kmer_dir, $basename . '.kmer');
        my $loc_file  = catfile($kmer_dir, $basename . '.loc');

        if (-e $kmer_file && -e $loc_file) {
            say "kmer/loc files exist";
            next FILE;
        }

        say "kmerizing";
        my $fa = Bio::SeqIO->new(-file => $file);
        open my $kmer_fh, '>', $kmer_file;
        open my $loc_fh,  '>', $loc_file;
        
        my $i = 0;
        while (my $seq = $fa->next_seq) {
            my $sequence = $seq->seq;
            my $num_kmers = length($sequence) + 1 - $mer_size;

            if ($num_kmers > 0) {
                for my $pos (0 .. $num_kmers - 1) {
                    say $kmer_fh join("\n",
                        '>' . $i++,
                        substr($sequence, $pos, $mer_size)
                    );
                }
            }

            print $loc_fh join("\t", $seq->id, $num_kmers), "\n"; 
        }
    }
}

# --------------------------------------------------
sub make_matrix {
    my $args       = shift;
    my $mode_dir   = catdir($args->{'out_dir'}, 'mode');
    my $matrix_dir = catdir($args->{'out_dir'}, 'matrix');

    unless (-d $mode_dir) {
        die "Bad mode dir ($mode_dir)";
    }

    unless (-d $matrix_dir) {
        make_path($matrix_dir);
    }

    my @files = File::Find::Rule->file()->in($mode_dir);

    unless (@files) {
        die "Cannot find any mode files in '$mode_dir'\n";
    }

    printf "Found %s mode files\n", scalar(@files);

    my %matrix;
    for my $file (@files) {
        open my $fh, '<', $file;
        chomp(my $n = <$fh>);
        close $fh;

        my $sample1 = basename(dirname($file));
        my $sample2 = basename($file);

        for ($sample1, $sample2) {
            $_ =~ s/\.\w+$//; # remove file extension
        }

        $matrix{ $sample1 }{ $sample2 } += $n || 0;
    }

#    my $num_samples = scalar(keys %matrix);
#    for my $sample (keys %matrix) {
#        if ($num_samples != values $matrix{ $sample }) {
#            say "Deleting $sample";
#            delete $matrix{ $sample };
#        }
#    }

    my @keys     = keys %matrix;
    my @all_keys = sort(uniq(@keys, map { keys %{ $matrix{ $_ } } } @keys));

    open my $fh, '>', catfile($matrix_dir, 'matrix.tab');

    say $fh join "\t", '', @all_keys;
    for my $sample1 (@all_keys) {
        my @vals = map { $matrix{ $sample1 }{ $_ } || 0 } @all_keys;

        say $fh join "\t", 
            $sample1, 
            #map { sprintf('%.2f', $_ > 0 ? log($_) : 0) } @vals,
            map { sprintf('%.2f', $_ > 0 ? log($_) : 1) } @vals,
        ;
    }
}

# --------------------------------------------------
sub pairwise_cmp {
    my $args          = shift;
    my $mode_min      = $args->{'mode_min'};
    my $jf_idx_dir    = catdir( $args->{'out_dir'}, 'jf' );
    my $kmer_dir      = catdir( $args->{'out_dir'}, 'kmer' );
    my $read_mode_dir = catdir( $args->{'out_dir'}, 'read_mode' );
    my $mode_dir      = catdir( $args->{'out_dir'}, 'mode' );
    my $tmp_dir       = catdir( $args->{'out_dir'}, 'tmp' );

    unless (-d $jf_idx_dir) {
        die "Bad Jellyfish index dir ($jf_idx_dir)";
    }

    unless (-d $kmer_dir) {
        die "Bad kmer dir ($kmer_dir)";
    }

    for my $dir ($tmp_dir, $mode_dir, $read_mode_dir) {
        make_path($dir) unless -d $dir;
    }

    my @jf_files       = File::Find::Rule->file()->in($jf_idx_dir);
    my $num_jf_files   = scalar(@jf_files);
    my @loc_files      = File::Find::Rule->file()->name('*.loc')->in($kmer_dir);
    my $num_loc_files  = scalar(@loc_files);
    my @kmer_files     = 
        File::Find::Rule->file()->name('*.kmer')->in($kmer_dir);
    my $num_kmer_files = scalar(@kmer_files);

    printf "Found %s Jellyfish files, %s kmer files\n", 
        $num_jf_files, $num_kmer_files;

    if ($num_jf_files < 1) {
        say "Not enough files to perform pairwise comparison.";
        return;
    }

    if ($num_jf_files != $num_kmer_files) {
        say "Number of Jellyfish files must equal kmer files.\n";
        return;
    }

    if ($num_loc_files != $num_kmer_files) {
        say "Number of location files must equal kmer files.\n";
        return;
    }

    my @combos;
    for my $pair (combine(2, @jf_files)) {
        my ($s1, $s2) = @$pair;
        push @combos, [$s1, $s2], [$s2, $s1];
    }

    printf "Will perform %s comparisons\n", scalar(@combos);

    my $combo_num = 0;
    COMBO:
    for my $pair (@combos) {
        my ($s1, $s2)       = @$pair;
        my $base_jf_file    = basename($s1, '.jf');
        my $base_kmer_file  = basename($s2, '.jf');
        my $kmer_file       = catfile($kmer_dir,  $base_kmer_file . '.kmer');
        my $loc_file        = catfile($kmer_dir, $base_kmer_file . '.loc');
        my $sample_mode_dir = catdir($mode_dir, $base_jf_file);
        my $sample_read_dir = catdir($read_mode_dir, $base_jf_file);
        my $mode_file       = catfile($sample_mode_dir, $base_kmer_file);
        my $read_mode_file  = catfile($sample_read_dir, $base_kmer_file);

        for my $dir ($sample_mode_dir, $sample_read_dir) {
            make_path($dir) unless -d $dir;
        }

        printf "%5d: %s -> %s", ++$combo_num, $base_kmer_file, $base_jf_file;

        if (-s $mode_file) {
            say " mode file exists";
            next COMBO;
        }
        else {
            say '';
        }

        my ($tmp_fh, $jf_query_out_file) = tempfile(DIR => $tmp_dir);
        close $tmp_fh;

        sys_exec('jellyfish', 'query', '-s', $kmer_file, 
            '-o', $jf_query_out_file, $s1);

        open my $loc_fh ,      '<', $loc_file;
        open my $mode_fh,      '>', $mode_file;
        open my $read_mode_fh, '>', $read_mode_file;
        open my $jf_fh,        '<', $jf_query_out_file;

        my $mode_count = 0;
        while (my $loc = <$loc_fh>) {
            chomp($loc);
            my ($read_id, $n_kmers) = split /\t/, $loc;

            my @counts;
            for my $val (take($n_kmers, $jf_fh)) {
                next if !$val;
                my ($kmer_seq, $count) = split /\s+/, $val;
                push @counts, $count if defined $count && $count =~ /^\d+$/;
            }

            my $mode = mode(@counts) // 0;
            if ($mode >= $mode_min) {
                print $read_mode_fh join("\t", $read_id, $mode), "\n";
                $mode_count++;
            }
        }

        say $mode_fh $mode_count;

        unlink $jf_query_out_file;
    }
}

# --------------------------------------------------
sub subset_files {
    my $args        = shift;
    my $max_seqs    = $args->{'max_seqs'};
    my $max_samples = $args->{'max_samples'};
    my @files       = File::Find::Rule->file()->in($args->{'in_dir'});

    printf "Found %s files in %s\n", scalar(@files), $args->{'in_dir'};

    if (scalar(@files) > $max_samples) {
        say "Subsetting to $max_samples files";
        @files = sample(-set => \@files, -sample_size => $max_samples);
    }

    $args->{'files'} = \@files;

    my $subset_dir = catdir($args->{'out_dir'}, 'subset');
    unless (-d $subset_dir) {
        make_path($subset_dir);
    }

    my $file_num = 0;
    FILE:
    for my $file (@files) {
        my $basename    = basename($file);
        my $subset_file = catfile($subset_dir, $basename);

        printf "%5d: %s, ", ++$file_num, $basename;

        if (-e $subset_file) {
            say "subset file exists";
            next FILE;
        }

        my ($tmp_fh, $tmp_filename) = tempfile();
        my $fa = Bio::SeqIO->new(-file => $file);
        my $count = 0;
        while (my $seq = $fa->next_seq) {
            $count++;
            say $tmp_fh $seq->id;
        }
        close $tmp_fh;

        print "$count seqs, ";

        if ($count < $max_seqs) {
            say "copying to subset file";
            copy($file, $subset_file); 
        }
        else {
            say "randomly sampling";

            my $random = File::RandomLine->new(
                $tmp_filename, 
                { algorithm => 'uniform' }
            ); 

            my %take;
            while (scalar(keys %take) < $max_seqs) {
                my $id = $random->next;
                $take{ $id }++;
            }

            my $in = Bio::SeqIO->new(-file => $file);
            my $out= Bio::SeqIO->new( 
                -format => 'Fasta', 
                -file => ">$subset_file"
            );

            while (my $seq = $in->next_seq) {
                $out->write_seq($seq) if exists $take{ $seq->id };
            }
        }

        unlink($tmp_filename);
    }

    return 1;
}

# --------------------------------------------------
sub mode {
    my @vals = @_ or return;
    my $mode = 0;

    if (scalar @vals == 1) {
        $mode = shift @vals;
    }
    else {
        my @distinct = uniq(@vals);

        if (scalar @distinct == 1) {
            $mode = shift @distinct;
        }
        else {
            my $stats = Statistics::Descriptive::Discrete->new;
            $stats->add_data(@vals);
            return $stats->mode();

#            if (my $mean = int($stats->mean())) {
#                my $two_stds = 2 * (int $stats->standard_deviation());
#                my $min      = $mean - $two_stds;
#                my $max      = $mean + $two_stds;
#
#                if (my @filtered = grep { $_ >= $min && $_ <= $max } @vals) {
#                    my $stats2 = Statistics::Descriptive::Discrete->new;
#                    $stats2->add_data(@filtered);
#                    $mode = int($stats2->mode());
#                }
#            }
#            else {
#                return 0;
#            }
        }
    }

    return $mode;
}

# ----------------------------------------------------
sub take {
    my ($n, $fh) = @_;

    my @return;
    for (my $i = 0; $i < $n; $i++) {
        my $line = <$fh>;
        last if !defined $line;
        chomp($line);
        push @return, $line;
    }

    @return;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

run-fizkin - run fizkin

=head1 SYNOPSIS

  run-fizkin -i input-dir -o output-dir

Required arguments:

  --input_dir   Input directory (FASTA)
  --output_dir  Output directory (FASTA)

Options (defaults in parentheses):
 
  --kmer_size    Kmer size (20)
  --mode_min     Minimum mode to take a sequence (1)
  --num_threads  Number of threads to use for Jellyfish (12)
  --hash_size    Size of hash for Jellyfish (100M)
  --max_seqs     Maximum number of sequences per input file (300,000)
  --max_samples  Maximum number of samples (15)
  --help         Show brief help and exit
  --man          Show full documentation

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

=head1 SEE ALSO

perl.

=head1 AUTHOR

Ken Youens-Clark E<lt>kyclark@email.arizona.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2015 kyclark

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
