#!/usr/bin/env perl
package mode::chip;

use Modern::Perl;
use MooseX::App::Command;
extends qw(mode);
use File::Path qw(remove_tree);
use File::Cat;
use Function;
use validate_options;
use input;

command_short_description q[Analysis for ChIP-seq];
command_long_description q[Analysis for ChIP-seq];
command_usage q[pRNASeqTools chip [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT1]=[file1]+[file2] ... --treatment [TREATMENT2]=[file1]+[file2] ... ];

option 'no-mapping' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Just perform the statistic analysis],
);
option 'mapping-only' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Do not perform the statistic analysis],
);
option 'foldchange' => (
  is => 'rw',
  isa => 'Num',
  default => '2',
  documentation => q[Threshold for differential peaks in fold enrichment.],
);
option 'qvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '0.05',
  documentation => q[Threshold for differential peaks in FDR.],
);
option 'igg' => (
  is => 'rw',
  isa => 'HashRef',
  documentation => q[IgG background.],
);
option 'style' => (
  is => 'rw',
  isa => 'Str',
  default => 'histone',
  documentation => q[Peak style, supported: factor, histone, tss.],
);
option 'control' => (
  is => 'rw',
  isa => 'HashRef',
  documentation => q[seq files seperated by plus for control],
);
option 'treatment' => (
  is => 'rw',
  isa => 'HashRef',
  required => 1,
  documentation => q[seq files seperated by plus for treatment, multiple treatments allowed],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $foldchange = $options{'foldchange'};
  my $qvalue = $options{'qvalue'};
  my $nomapping = $options{'no-mapping'};
  my $mappingonly = $options{'mapping-only'};
  my $input = $options{'control'};
  my $ip = $options{'treatment'};
  my $igg = $options{'igg'};
  my $style = $options{'style'};

  my (@tags, @files);
  my ($input_ref, $inputf_ref, $inputp_ref, $ip_ref, $ipf_ref, $ipp_ref, $igg_ref, $iggf_ref, $iggp_ref);
  my $homer_input = "";
  my $homer_bg = "";
  if(defined $input){
    ($input_ref, $inputf_ref, $inputp_ref) = input->run($input);
    push @tags, @$input_ref;
    push @files, @$inputf_ref;
    $homer_input = "-i ".join " ", @$input_ref;
  }
  ($ip_ref, $ipf_ref, $ipp_ref) = input->run($ip);
  push @tags, @$ip_ref;
  push @files, @$ipf_ref;
  my @ipp = @$ipp_ref;
  my $homer_ip = "-t ".join " ", @$ip_ref;
  if (defined $igg){
    ($igg_ref, $iggf_ref, $iggp_ref) = input->run($igg);
    push @tags, @$igg_ref;
    push @files, @$iggf_ref;
    $homer_bg = "-b ".join " ", @$igg_ref;
  }

  if(!$nomapping){
    print $main::tee "\nBuilding index...\n";
    symlink $prefix."/reference/".$genome."_chr_all.fasta", $genome."_chr_all.fasta";
    system ("bowtie2-build -q ".$genome."_chr_all.fasta ".$genome."_chr_all");
    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];
  		print $main::tee "\nMapping $tag...\n";
  		if($file !~ /,/){
  			my @files = Function->SRR($file, $thread);
        if($#files == 0){
        	Function->unzip($files[0], $tag);
        	if(defined $adaptor){
            print $main::tee "\nTrimming...\n";
						system ("cutadapt -j ".$thread." -m 50 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
          	rename $tag."_trimmed.fastq", $tag.".fastq";
        	}
    			system ("bowtie2 -p ".$thread." -x ".$genome."_chr_all -U ".$tag.".fastq -S ".$tag.".sam 2>&1");
    			unlink ($tag.".fastq");
        }else{
          Function->unzip($files[0], $tag."_R1");
          Function->unzip($files[1], $tag."_R2");
          if(defined $adaptor){
            system ("cutadapt -j ".$thread." -m 50 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
            rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
            rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
          }
          system ("bowtie2 -p ".$thread." -x ".$genome."_chr_all -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq -S ".$tag.".sam 2>&1");
          unlink ($tag."_R1.fastq", $tag."_R2.fastq");
        }
  		}else{
       	my ($file1, $file2) = split /,/, $file;
				Function->unzip($file1, $tag."_R1");
				Function->unzip($file2, $tag."_R2");
  		  if(defined $adaptor){
          system ("cutadapt -j ".$thread." -m 50 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
          rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
        }
        system ("bowtie2 -p ".$thread." -x ".$genome."_chr_all -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq -S ".$tag.".sam 2>&1");
  			unlink ($tag."_R1.fastq", $tag."_R2.fastq");
  		}
  		system ("samtools view -Sb -q 10 --threads ".$thread." ".$tag.".sam > ".$tag.".bam");
      system ("samtools sort -o ".$tag.".sorted.bam ".$tag.".bam");
      system ("samtools index ".$tag.".sorted.bam");
      unlink $tag.".sam", $tag.".bam";
      system ("bamCoverage -b ".$tag.".sorted.bam -bs 5 -p ".$thread." --ignoreDuplicates --normalizeUsing CPM -o ".$tag.".bw");
  		print $main::tee "\nMapping completed!\n";
      system ("makeTagDirectory ".$tag." ".$tag.".sorted.bam -tbp 1 -genome ".$genome."_chr_all.fasta 2>&1");
    }
    unlink (glob ($genome."_chr_all*"), "igv.log");
    if(!$mappingonly && defined $input){
      print $main::tee "\nFinding peaks...\nFold Enrichment\t$foldchange\tQ Value\t$qvalue\n";
      system ("getDifferentialPeaksReplicates.pl -genome ".$prefix."/reference/".$genome."_chr_all.fasta -gff ".$prefix."/reference/".$genome."_genes.gff -f ".$foldchange." -q ".$qvalue." -tbp 1 -style ".$style." -fdr 0.01 -F 2 -L 2 ".$homer_input." ".$homer_bg." ".$homer_ip." > ".$ipp[0].".".$style.".peaks.txt");
    }
	}else{
    foreach my $pre (@tags){
      symlink "../".$pre, $pre or die $!;
    }
    print $main::tee "\nFinding peaks...\nFold Enrichment\t$foldchange\tQ Value\t$qvalue\n";
    system ("getDifferentialPeaksReplicates.pl -genome ".$prefix."/reference/".$genome."_chr_all.fasta -gff ".$prefix."/reference/".$genome."_genes.gff -f ".$foldchange." -q ".$qvalue." -tbp 1 -style ".$style." -fdr 0.01 -F 2 -L 2 ".$homer_input." ".$homer_bg." ".$homer_ip." > ".$ipp[0].".".$style.".peaks.txt");
    foreach my $pre (@tags){
      remove_tree $pre;
    }
  }
}

1;

__END__
