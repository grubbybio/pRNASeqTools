#!/usr/bin/env perl
package mode::atac;

use Modern::Perl;
use MooseX::App::Command;
extends qw(mode);
use File::Path qw(remove_tree);
use File::Cat;
use Function;
use validate_options;
use input;

command_short_description q[Analysis for ATAC-seq];
command_long_description q[Analysis for ATAC-seq];
command_usage q[pRNASeqTools atac [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT1]=[file1]+[file2] ... --treatment [TREATMENT2]=[file1]+[file2] ... ];

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
option 'auc' => (
  is => 'rw',
  isa => 'Num',
  default => '20',
  documentation => q[Threshold for peaks in area under curve.],
);
option 'qvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '1',
  documentation => q[Threshold for differential peaks in FDR.],
);
option 'pvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '0.01',
  documentation => q[Threshold for differential peaks in FDR.],
);
option 'control' => (
  is => 'rw',
  isa => 'HashRef',
  documentation => q[seq files seperated by plus for control, either input or IgG],
);
option 'treatment' => (
  is => 'rw',
  isa => 'HashRef',
  required => 1,
  documentation => q[seq files seperated by plus for IP],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $auc = $options{'auc'};
  my $qvalue = $options{'qvalue'};
  my $pvalue = $options{'pvalue'};
  my $nomapping = $options{'no-mapping'};
  my $mappingonly = $options{'mapping-only'};
  my $input = $options{'control'};
  my $ip = $options{'treatment'};

  my (@tags, @files);
  my ($input_ref, $inputf_ref, $inputp_ref, $ip_ref, $ipf_ref, $ipp_ref);
  my $genrich_input = "";
  my $genrich_ip = "";
  if(defined $input){
    ($input_ref, $inputf_ref, $inputp_ref) = input->run($input);
    push @tags, @$input_ref;
    push @files, @$inputf_ref;
    my @inputp = @$inputp_ref;
    $genrich_input = "-c ".join ",", map {$_.".sorted.name.bam"} @$input_ref;
  }
  ($ip_ref, $ipf_ref, $ipp_ref) = input->run($ip);
  push @tags, @$ip_ref;
  push @files, @$ipf_ref;
  my @ipp = @$ipp_ref;
  $genrich_ip = "-t ".join ",", map {$_.".sorted.name.bam"} @$ip_ref;

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
      unlink $tag.".sam";
      system ("samtools sort -n -o ".$tag.".sorted.name.bam ".$tag.".bam");
      system ("samtools fixmate -m ".$tag.".sorted.name.bam ".$tag.".fixmate.bam");
      system ("samtools sort -o ".$tag.".sorted.bam ".$tag.".fixmate.bam");
      system ("samtools markdup -r ".$tag.".sorted.bam ".$tag.".sorted.dedup.bam");
      system ("samtools index ".$tag.".sorted.dedup.bam");
      unlink $tag.".bam", $tag.".fixmate.bam", $tag.".sorted.bam";
      system ("bamCoverage -b ".$tag.".sorted.dedup.bam --skipNAs -bs 5 -p ".$thread." --ignoreDuplicates --normalizeUsing CPM -o ".$tag.".bw");
  		print $main::tee "\nMapping completed!\n";
      system ("Genrich -r -v -g 200 -p 0.001 -t ".$tag.".sorted.name.bam -o ".$tag.".narrowPeak.txt");
    }
    unlink (glob ($genome."_chr_all*"), "igv.log");
    if(!$mappingonly && defined $input){
      my $command = "Genrich -j -r -v -a 20 -e chrC,chrM ".$genrich_ip." ".$genrich_input." -o ".$ipp[0].".gain.narrowPeak.txt";
      if($qvalue < 1){
        print $main::tee "\nFinding peaks...\nAUC\t$auc\tQ Value\t$qvalue\n";
        $command .= " -q $qvalue";
      }else{
        print $main::tee "\nFinding peaks...\nAUC\t$auc\tP Value\t$pvalue\n";
        $command .= " -p $pvalue";
      }
      system $command." 2>&1";
    }
	}else{
    foreach my $pre (@tags){
      symlink "../".$pre.".sorted.name.bam", $pre.".sorted.name.bam" or die $!;
    }
    my $command = "Genrich -j -r -v -a 20 ".$genrich_ip." ".$genrich_input." -o ".$ipp[0].".gain.narrowPeak.txt";
    if($qvalue < 1){
      print $main::tee "\nFinding peaks...\nAUC\t$auc\tQ Value\t$qvalue\n";
      $command .= " -q $qvalue";
      print $main::tee "$command\n";
    }else{
      print $main::tee "\nFinding peaks...\nAUC\t$auc\tP Value\t$pvalue\n";
      $command .= " -p $pvalue";
    }
    system $command." 2>&1";
    foreach my $pre (@tags){
      remove_tree $pre;
    }
  }
}

1;

__END__
