#!/usr/bin/env perl
package mode::wgbs;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use Function;
use File::Path qw/remove_tree/;
use Ref;
use validate_options;
use input;
use mode::wgbs;

command_short_description q[Analysis for WGBS sequencing];
command_long_description q[Analysis for WGBS sequencing];
command_usage q[pRNASeqTools wgbs [OPTIONS] --control [CONTROL]=[file1]+[file2] ... ];

option 'no-mapping' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Just perform the ribo-seq analysis],
);
option 'mapping-only' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Do not perform the ribo-seq analysis and only align the reads to the transcriptome],
);
option 'binsize' => (
  is => 'rw',
  isa => 'Num',
  default => 100,
  documentation => q[Windows size to slice the genome.],
);
option 'minC' => (
  is => 'rw',
  isa => 'Num',
  default => 4,
  documentation => q[Minimum reads per cytosine allowed.],
);

sub run {
  my ($self) = @_;
  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $nomapping = $options{'no-mapping'};
  my $mappingonly = $options{'mapping-only'};
  my $control = $options{'control'};
  my $treatment = $options{'treatment'};
  my $ns = $options{'ns'};
  my $h = $options{'h'};
  my $binsize = $options{'binsize'};
  my $minC = $options{'minC'};

  my ($tags_ref, $files_ref, $par_ref) = input->run($control);
  my @tags = @$tags_ref;
  my @files = @$files_ref;
  my @par = @$par_ref;
  if (defined $treatment){
    ($tags_ref, $files_ref, $par_ref) = input->run($treatment);
    push @tags, @$tags_ref;
    push @files, @$files_ref;
    push @par, @$par_ref;
  }
  my $par = join " ", @par;

  if(!$nomapping){
    print $main::tee "\nBuilding indices...\n";
    symlink $prefix."/reference/".$genome."_chr_all.fasta", $genome.".fasta";
    system ("bismark_genome_preparation . 2>&1");

    for(my $i=0;$i<=$#tags;$i++){
      my $file = $files[$i];
      my $tag = $tags[$i];

      print $main::tee "\nMapping $tag...\n";

      if($file !~ /,/){
        my @files = Function->SRR($file, $thread);
        if($#files == 0){
          Function->unzip($files[0], $tag);
          if(defined $adaptor){

            print $main::tee "\nStart trimming...\r";

            system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
            rename $tag."_trimmed.fastq", $tag.".fastq";
          }
          system ("bismark -p ".$thread." -N 1 . ".$tag.".fastq 2>&1");
          system ("deduplicate_bismark -s --bam ".$tag."_bismark_bt2.bam 2>&1");
          rename $tag."_bismark_bt2.deduplicated.bam", $tag.".bam";
          system ("bismark_methylation_extractor --parallel ".$thread." -s --bedGraph --cutoff 4 --cytosine_report --CX --genome_folder . ".$tag.".bam 2>&1");
          unlink ($tag.".fastq", $tag."_bismark_bt2.bam");
        }else{
          Function->unzip($files[0], $tag."_R1");
          Function->unzip($files[1], $tag."_R2");
          if(defined $adaptor){
            system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
            rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
            rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
          }
          system ("bismark -p ".$thread." -N 1 . -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq 2>&1");
          system ("deduplicate_bismark -p --bam ".$tag."_R1_bismark_bt2_pe.bam 2>&1");
          rename $tag."_R1_bismark_bt2_pe.deduplicated.bam", $tag.".bam";
          system ("bismark_methylation_extractor --parallel ".$thread." -p --bedGraph --cutoff 4 --cytosine_report --CX --genome_folder . ".$tag.".bam", $tag.".bam 2>&1");
          unlink ($tag."_R1.fastq", $tag."_R2.fastq", $tag."_R1_bismark_bt2_pe.bam")
        }
      }else{
        my ($file1, $file2) = split /,/, $file;
        Function->unzip($file1, $tag."_R1");
        Function->unzip($file2, $tag."_R2");
        if(defined $adaptor){
          system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
          rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
        }
        system ("bismark -p ".$thread." -N 1 . -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq 2>&1");
        system ("deduplicate_bismark -p --bam ".$tag."_R1_bismark_bt2_pe.bam 2>&1");
        rename $tag."_R1_bismark_bt2_pe.deduplicated.bam", $tag.".bam";
        system ("bismark_methylation_extractor --parallel ".$thread." -p --bedGraph --cutoff 4 --cytosine_report --CX --genome_folder . ".$tag.".bam 2>&1");
        unlink ($tag."_R1.fastq", $tag."_R2.fastq", $tag."_R1_bismark_bt2_pe.bam");
        system ("awk '{OFS=\"\\t\";if(\$4+\$5>0){if(\$6==\"CG\"){print \$1,\$2,\$2+1,\$4/(\$4+\$5) > \"".$tag.".CG.bedgraph\"}; if(\$6==\"CHG\"){print \$1,\$2,\$2+1,\$4/(\$4+\$5) > \"".$tag.".CHG.bedgraph\"}; if(\$6==\"CHH\"){print \$1,\$2,\$2+1,\$4/(\$4+\$5) > \"".$tag.".CHH.bedgraph\"}}}' ".$tag.".CX_report.txt");
        system ("sort -k1,1 -k2,2n ".$tag.".CX_report.txt > tmp");
        system ("bgzip -c tmp > ".$tag.".CX_report.txt.gz");
        system ("tabix -C -p vcf ".$tag.".CX_report.txt.gz");
      }
      print $main::tee "\nAlignment finished...\n";
      mode::wgbs->bin($tag, $binsize, $minC);
      unlink glob ("C*_O?_".$tag.".txt"), "tmp";
    }

    if(!$mappingonly and $#par > 1){
      print $main::tee "\nPerforming DMRcaller...\n";
      system ("Rscript --vanilla ".$prefix."/scripts/DMRcaller.R ".$thread." ".$par);
    }
    remove_tree "Bisulfite_Genome";
    unlink ($genome.".fasta", glob ("*.ebwt"));

  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".CX_report.txt.gz", $pre.".CX_report.txt.gz" or die $!;
    }
    print $main::tee "\nPerforming DMRcaller...\n";
    system ("Rscript --vanilla ".$prefix."/scripts/DMRcaller.R ".$thread." ".$par);
    unlink glob ("*.CX_report.txt.gz");
  }
}

sub bin {
  use POSIX;
  my ($self, $tag, $binsize, $minC) = @_;
  my @context = qw/CG CHG CHH/;
  my @filehandle = (*CG, *CHG, *CHH);
  my (%ccounttotal, %ctcounttotal, %cavg, %num, %covtotal, %sum, %totalbin, %passedbin, $chr, $flag, %ccount, %ctcount);
  foreach my $type (@context){
    $ccounttotal{$type} = 0;
    $ctcounttotal{$type} = 0;
    $cavg{$type} = 0;
    $num{$type} = 0;
    $covtotal{$type} = 0;
    $sum{$type} = 0;
    $totalbin{$type} = 0;
    $passedbin{$type} = 0;
  }
  my $currchr = "";

  foreach my $cont (0 .. $#context){
    open $filehandle[$cont], ">".$tag.".bin.".$binsize.".".$context[$cont].".txt";
    $filehandle[$cont]->print("bin\tccount\tctcount\tno.cytosin\tno.qualified.coverage\n");
  }

  open MAP, "<".$tag.".CX_report.txt" or die $!;
  while (my $row = <MAP>) {
    chomp $row;
    my @data = split /\t/, $row;
    if($. == 1){
      $chr =  $data[0];
      $flag = $binsize;

      print $main::tee "Working on $chr...\n";

    }
    my $currchr = $data[0];
    if($currchr ne $chr){
      foreach my $cont (0 .. $#context){
        my $avgMeth = 0;
        if ($num{$context[$cont]} > 0) {
          $avgMeth = $cavg{$context[$cont]} / $num{$context[$cont]};
          my $totalReads = $ccounttotal{$context[$cont]} + $ctcounttotal{$context[$cont]};
          $ccounttotal{$context[$cont]} = floor($totalReads * $avgMeth + 0.5);
          $ctcounttotal{$context[$cont]} = floor($totalReads * (1 - $avgMeth) + 0.5);
        }
        my $bin = $flag / 100;
        $filehandle[$cont]->print("$chr\_$bin\t$ccounttotal{$context[$cont]}\t$ctcounttotal{$context[$cont]}\t$num{$context[$cont]}\t$covtotal{$context[$cont]}\n");
      }
      $flag = $binsize;
      $chr = $currchr;

      print $main::tee "Working on $chr...\n";

      next;
    }
    my $pos = $data[1];
    my $strand = $data[2];
    my $type = $data[5];
    $ccount{$type} = $data[3];
    $ctcount{$type} = $data[4];
    $sum{$type} ++;

    while ($pos >= $flag){
      foreach my $cont (0 .. $#context){
        my $avgMeth = 0;
        if ($num{$context[$cont]} > 0) {
          $avgMeth = $cavg{$context[$cont]} / $num{$context[$cont]};
          my $totalReads = $ccounttotal{$context[$cont]} + $ctcounttotal{$context[$cont]};
          $ccounttotal{$context[$cont]} = floor($totalReads * $avgMeth + 0.5);
          $ctcounttotal{$context[$cont]} = floor($totalReads * (1 - $avgMeth) + 0.5);
        }
        my $bin = $flag / 100;
        $filehandle[$cont]->print("$chr\_$bin\t$ccounttotal{$context[$cont]}\t$ctcounttotal{$context[$cont]}\t$num{$context[$cont]}\t$covtotal{$context[$cont]}\n");
        $totalbin{$context[$cont]} ++;
        $passedbin{$context[$cont]} ++ if($covtotal{$context[$cont]} >= 4);

        $ccounttotal{$context[$cont]} = 0;
        $ctcounttotal{$context[$cont]} = 0;
        $cavg{$context[$cont]} = 0;
        $num{$context[$cont]} = 0;
        $covtotal{$context[$cont]} = 0;
      }
      $flag += $binsize;
    }

    next if($ccount{$type} < 1 && $ctcount{$type} < 1);
    next if($ccount{$type} + $ctcount{$type} < $minC);
    $ccounttotal{$type}  +=  $ccount{$type};
    $ctcounttotal{$type} +=  $ctcount{$type};
    $num{$type} ++;
    my $curAvg = $ccount{$type} / ($ccount{$type} + $ctcount{$type});
    $cavg{$type} += $curAvg;
    $covtotal{$type} ++ if ($ctcount{$type} + $ccount{$type}>= 4);
  }
  foreach my $cont (0 .. $#context){
    my $avgMeth = 0;
    if ($num{$context[$cont]} > 0) {
      $avgMeth = $cavg{$context[$cont]}/$num{$context[$cont]};
      my $totalReads = $ccounttotal{$context[$cont]}+$ctcounttotal{$context[$cont]};
      $ccounttotal{$context[$cont]} = floor($totalReads * $avgMeth + 0.5);
      $ctcounttotal{$context[$cont]} = floor($totalReads * (1 - $avgMeth) + 0.5);
    }
    my $bin = $flag / 100;
    $filehandle[$cont]->print("$chr\_$bin\t$ccounttotal{$context[$cont]}\t$ctcounttotal{$context[$cont]}\t$num{$context[$cont]}\t$covtotal{$context[$cont]}\n");
  }
  foreach my $reftype (@context){

    print $main::tee "$tag\t$reftype\t$totalbin{$reftype}\t$passedbin{$reftype}\n";

  }

  close MAP;
  close CG;
  close CHG;
  close CHH;
}

1;
__END__
