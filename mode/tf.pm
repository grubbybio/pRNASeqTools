#!/usr/bin/env perl
package mode::tf;

use Modern::Perl;
use MooseX::App::Command;
extends qw(mode);
use File::Path qw(remove_tree make_path);
use File::Cat;
use Function;
use Ref;
use validate_options;
use input;

command_short_description q[Analysis for two-factor DE];
command_long_description q[Analysis for two-factor DE identification following modes "mrna" and "srna"];
command_usage q[pRNASeqTools tf [OPTIONS] --control [CONTROL]=[DESIGN1],[N1],[DESIGN2],[N2]... --treatment [TREATMENT]=[DESIGN1],[N1],[DESIGN2],[N2]...];

option 'foldchange' => (
  is => 'rw',
  isa => 'Num',
  default => '1.5',
  documentation => q[Threshold for DEG in fold change.],
);
option 'pvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '0.05',
  documentation => q[Threshold for DEG in P value.],
);
option 'mode' => (
  is => 'rw',
  isa => 'Str',
  default => 'mrna',
  documentation => q[Data from transcriptome or sRNAome, supported values: mrna, srna.],
);
option 'norm' => (
  is => 'rw',
  isa => 'Str',
  default => 'rRNA,total',
  documentation => q[Method for normalization, seperated by comma. Allowed: rRNA, total, spike_in.],
);
option 'binsize' => (
  is => 'rw',
  isa => 'Num',
  default => '100',
  documentation => q[Windows size to slice the genome.],
);
option 'DESeq2Norm' => (
  is => 'rw',
  isa => 'Str',
  default => 'DESeq2',
  documentation => q[Normalization method for DESeq2, "DESeq2" and "RPM" are allowed],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $prefix = $options{'prefix'};
  my $foldchange = $options{'foldchange'};
  my $pvalue = $options{'pvalue'};
  my $control = $options{'control'};
  my $treatment = $options{'treatment'};
  my $mode = $options{'mode'};
  my $norm = $options{'norm'};
  my @norms = split /,/, $norm;
  my $binsize = $options{'binsize'};
  my $genome = $options{'genome'};
  my $deseq2norm = $options{'DESeq2Norm'};

  print $main::tee "Two-factor comparison between ".$control." and ".$treatment."\nFold change = ".$foldchange." P value = ".$pvalue."\n";

  my ($tag, $rep, @tags);
  my $n = 0;
  my ($control_genotype) = keys %$control;
  my ($treatment_genotype) = keys %$treatment;
  my @controls = split /,/, $$control{$control_genotype};
  my @treatments = split /,/, $$treatment{$treatment_genotype};

  die "Please provide paired data!" if ($#controls != $#treatments);

  my $pars = join " ", ($control_genotype, @controls, $treatment_genotype, @treatments);
  for(@controls){
    ($tag, $rep) = splice(@controls, 0, 2);
    for($n=1;$n<=$rep;$n++){
      push @tags, $control_genotype."_".$tag."_".$n;
    }
  }
  for(@treatments){
    ($tag, $rep) = splice(@treatments, 0, 2);
    for($n=1;$n<=$rep;$n++){
      push @tags, $treatment_genotype."_".$tag."_".$n;
    }
  }
  if($mode eq 'mrna'){
    foreach my $pre (@tags){
      symlink "../".$pre.".txt", $pre.".txt" or die $!;
    }
    system("Rscript --vanilla ".$prefix."/scripts/tf_mrna.R ".$deseq2norm." ".$pvalue." ".$foldchange." ".$pars);
    unlink glob "*_?.txt";
  }elsif($mode eq 'srna'){
    foreach my $pre (@tags){
      symlink "../".$pre.".nf", $pre.".nf" or die $!;
      opendir my $originalFolder, "../" or die $!;
      my @fileList = readdir $originalFolder;
      foreach my $file (@fileList){
        symlink "../".$file, $file if($file =~ /^$pre.+count$/ and $file !~ /norm/);
      }
      closedir $originalFolder;
    }
    foreach my $mnorm (@norms){
      system("Rscript --vanilla ".$prefix."/scripts/tf_srna.R ".$mnorm." ".$pvalue." ".$foldchange." ".$pars);
      opendir my $dir, "." or die $!;
      my @dirdd = grep {/csv$/} readdir $dir;
      my @dird = grep {/$mnorm\.bin/} @dirdd;
      my @dir = grep {/hyper|hypo/} @dird;
      closedir $dir;
      foreach my $hcsv (@dir){
        open CSV, "$hcsv" or die $!;
        (my $bg = $hcsv) =~ s/csv$/bedgraph/;
        open TMP, ">$bg" or die $!;
        while (my $ff = <CSV>){
          chomp $ff;
          my @row = split /,/, $ff;
          if($row[0]=~/(\w+)_(\d+)/){
            my $chr = $1;
            my $start = $2 * $binsize;
            my $end = $start + $binsize - 1;
            print TMP "$chr\t$start\t$end\t$row[2]\n";
          }
        }
        close TMP;
        close CSV;
      }
      my %ann = Ref->ann($prefix, $genome, $binsize);
    	foreach my $csv (@dird){
    		open CSV, "$csv" or die $!;
    		open TMP, ">tmp4" or die $!;
    		while(my $bb = <CSV>){
    			chomp $bb;
    			my @row = split /,/, $bb;
    			$row[0] =~ s/"//g;
    			if(exists $ann{$row[0]}){
    				print TMP "$bb,$ann{$row[0]}\n";
    			}else{
    				print TMP "$bb,Intergenic\n";
    			}
    		}
    		close CSV;
    		close TMP;
    		rename "tmp4", $csv;
    	}
      system("Rscript --vanilla ".$prefix."/scripts/tf_mirna.R ".$mnorm." ".$pvalue." ".$foldchange." ".$pars);
    }
    unlink glob ("*.count"), glob ("*.nf");
  }
}

1;
