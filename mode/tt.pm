#!/usr/bin/env perl
package mode::tt;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use Function;
use Ref;
use File::Path qw/remove_tree/;
use validate_options;
use input;

command_short_description q[Analysis for trunction and tailing of miRNAs];
command_long_description q[Analysis for trunction and tailing of miRNAs];
command_usage q[pRNASeqTools tt [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT]=[file1]+[file2] ... ];

option 'mmap' => (
  is => 'rw',
  isa => 'Str',
  default => 'u',
  documentation => q[],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $mmap = $options{'mmap'};
  my $control = $options{'control'};
  my $treatment = $options{'treatment'};
  my %fas = Ref->fas($prefix, $genome);
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

  for(my $i=0;$i<=$#tags;$i++){
    my $file = $files[$i];
    my $tag = $tags[$i];

    print $main::tee "\nMapping $tag...\n";

    $file = Function->SRR($file, $thread);
    Function->unzip($file, $tag);
    if(defined $adaptor){

      print $main::tee "\nStart trimming...\n";

      system ("cutadapt -j ".$thread." -m 14 -M 42 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
      rename $tag."_trimmed.fastq", $tag.".fastq";
    }
    system ("bowtie -v 0 -p ".$thread." -t --un ".$tag.".unmapped_0.fastq --al ".$tag.".mapped.fastq ".$prefix."/reference/".$genome."_chr_all ".$tag.".fastq ".$tag.".out 2>&1");
    my ($seq, $plus, $qua, $id);
    for(my $p = 1; $p <= 8; $p ++){
      my $j = $p - 1;
      my $ss = "";
      open FQ, "<$tag.unmapped_$j.fastq" or die $!;
      open OUT, ">tmp.fastq";
      my $n = 1;
      while(my $aa = <FQ> ){
        chomp $aa;
        if($n == 4){
          $qua = substr($aa, 0, -1);
        if(length($seq) >= 14){
          print OUT "$id\_$ss\n$seq\n$plus\n$qua\n";
        }
        $n = 1;
        next;
        }elsif($n == 3){
          $plus = $aa;
          $n = 4;
        }elsif($n == 2){
          $seq = substr($aa, 0, -1);
          $ss .= lc(substr($aa, -1));
          $n = 3;
        }elsif($n == 1){
          my @row = split /\s/, $aa;
          ($id, $ss) = split /_/, $row[0];
          $n = 2;
        }
      }
      close OUT;
      close FQ;
      unlink ($tag.".unmapped_$j.fastq");
      system ("bowtie -v 0 -p ".$thread." -t --un ".$tag.".unmapped_".$p.".fastq --al ".$tag.".mapped_".$p.".fastq ".$prefix."/reference/".$genome."_chr_all tmp.fastq ".$tag.".out 2>&1");
    }
    system ("cat $tag.mapped*.fastq > ".$tag."_edited.fastq");
    unlink (glob ($tag."*mapped*"), "tmp.fastq");
    system ("ShortStack --outdir ".$tag."tmp --align_only --bowtie_m 1000 --ranmax 50 --mmap ".$mmap." --mismatches 0 --bowtie_cores ".$thread." --nohp --readfile ".$tag."_edited.fastq --genomefile ".$prefix."/reference/".$genome."_chr_all.fasta 2>&1");
    system ("samtools view -h ".$tag."tmp/$tag\_edited.bam > $tag");
    system ("awk '{if(\$10!=\"*\" && \$3!=\"*\") print > (FILENAME\".edited.sam\")}' ".$tag);
    system ("samtools view -Sb ".$tag.".edited.sam > ".$tag.".edited.bam");
    remove_tree $tag."tmp";
    unlink ($tag."_edited.fastq", $tag.".fastq", $tag, $tag.".edited.sam");

    print $main::tee "\nAlignment Completed!\n";

    my %mir = Ref->mir($prefix, $genome);
    system ("bedtools intersect -wo -s -a ".$prefix."/reference/".$genome."_miRNA_miRNA_star.gff -b ".$tag.".edited.bam > ".$tag.".out");
    open OV, "<$tag.out" or die $!;
    my ($name, $tail);
    while (my $cc = <OV>){
      chomp $cc;
      my @data = split /\t/, $cc;
      my @row = split /_/, $data[12];
      if($row[1]){
        $tail = reverse $row[1];
      }else{
        $tail = "";
      }
      $mir{$data[8]}{"read"}{$data[10]."_".$data[11]."_".$tail}{"count"} ++;
    }
    close OV;

    open OUTPUT1, ">$tag.seq.out" or die $!;
    my ($refseq, $miseq, $flank1, $flank2, %out);
    my $flank = 25;
    foreach my $mm (sort keys %mir){
      my $lim1 = $mir{$mm}{"start"} - $flank;
      my $lim2 = $mir{$mm}{"end"} + $flank;
      my $length = $mir{$mm}{"end"} - $mir{$mm}{"start"} + 1;
      $refseq = substr($fas{$mir{$mm}{"chromosome"}}, $lim1, $lim2-$lim1+1);
      $miseq = substr($fas{$mir{$mm}{"chromosome"}}, $mir{$mm}{"start"}-1, $mir{$mm}{"end"}-$mir{$mm}{"start"}+1);
      if($mir{$mm}{"strand"} eq "+"){
        print OUTPUT1 "$refseq\t$mm\n"."."x($flank-1).$miseq."."x($flank+1)."\tREF\n";
      }else{
        $refseq = Function->revcomp($refseq);
        $miseq = Function->revcomp($miseq);
        print OUTPUT1 "$refseq\t$mm\n"."."x($flank+1).$miseq."."x($flank-1)."\tREF\n";
      }
      foreach my $read (sort keys %{$mir{$mm}{"read"}}){
        my @features = split /_/, $read;
        $features[2] = "" if(!$features[2]);
        if($mir{$mm}{"strand"} eq "+"){
          $seq = substr($fas{$mir{$mm}{"chromosome"}}, $features[0], $features[1]-$features[0]);
          $flank1 = $features[0] - $lim1;
          $flank2 = $lim2 - $features[1] - length($features[2]);
          if($mir{$mm}{"start"} == $features[0]+1){
            if($length-$features[1]+$features[0]-1 < 0){
              $out{$mm}{$features[1]-$features[0]-$length+length($features[2])}{0} += $mir{$mm}{"read"}{$read}{"count"};
            }else{
              $out{$mm}{length($features[2])}{$length-$features[1]+$features[0]} += $mir{$mm}{"read"}{$read}{"count"};
            }
          }
        }else{
          $seq = substr($fas{$mir{$mm}{"chromosome"}}, $features[0], $features[1]-$features[0]);
          $seq = Function->revcomp($seq);
          $flank1 = $lim2 - $features[1] + 1;
          $flank2 = $features[0] - $lim1 - length($features[2]) - 1;
          if($mir{$mm}{"end"} == $features[1]){
            if($length-$features[1]+$features[0]-1 < 0){
              $out{$mm}{$features[1]-$features[0]-$length+length($features[2])}{0} += $mir{$mm}{"read"}{$read}{"count"};
            }else{
              $out{$mm}{length($features[2])}{$length-$features[1]+$features[0]} += $mir{$mm}{"read"}{$read}{"count"};
            }
          }
        }
        print OUTPUT1 "."x$flank1.$seq.$features[2]."."x($flank2+1)."\t$mir{$mm}{read}{$read}{count}\n";
      }
    }
    close OUTPUT1;
    open OUTPUT2, ">$tag.out" or die $!;
    foreach my $mm (sort keys %mir){
      for my $tailing (0 .. 8){
        for my $truncation (0 .. 8){
          if(exists $out{$mm}{$tailing}{$truncation}){
            print OUTPUT2 "$mm\t$tailing\t$truncation\t$out{$mm}{$tailing}{$truncation}\n";
          }else{
            print OUTPUT2 "$mm\t$tailing\t$truncation\t0\n";
          }
        }
      }
    }
    close OUTPUT2;

    open GFF, $prefix."/reference/".$genome."_MIR.gff" or die $!;
    my %miR =();
    while (my $bb = <GFF>){
      chomp $bb;
      my @row = split /\t/, $bb;
      $miR{$row[8]}{"chromosome"} = $row[0];
      $miR{$row[8]}{"strand"} = $row[6];
      $miR{$row[8]}{"start"} = $row[3];
      $miR{$row[8]}{"end"} = $row[4];
    }
    close GFF;
    system ("bedtools intersect -wo -s -a ".$prefix."/reference/".$genome."_MIR.gff -b ".$tag.".edited.bam > ".$tag.".out2");
    open OV, "<$tag.out2" or die $!;
    while (my $cc = <OV>){
      chomp $cc;
      my @data = split /\t/, $cc;
      my @row = split /_/, $data[12];
      if($row[1]){
        $tail = reverse $row[1];
      }else{
        $tail = "";
      }
      $miR{$data[8]}{"read"}{$data[10]."_".$data[11]."_".$tail}{"count"} ++;
    }
    close OV;

    open OUTPUT3, ">$tag.seq.out2" or die $!;
    my %out2 = ();
    foreach my $mm (sort keys %miR){
      my $lim1 = $miR{$mm}{"start"} - $flank;
      my $lim2 = $miR{$mm}{"end"} + $flank;
      my $length = $miR{$mm}{"end"} - $miR{$mm}{"start"} + 1;
      my @loc;
      my $miseq2 = "";
      $refseq = substr($fas{$miR{$mm}{"chromosome"}}, $lim1, $lim2-$lim1+1);
      foreach my $mip (sort keys %mir){
        if($mip =~ /$mm/i){
          push @loc, $mir{$mip}{"start"};
          push @loc, $mir{$mip}{"end"};
        }
      }
      @loc = sort {$a <=> $b} @loc;
      my $miseq1 = substr($fas{$miR{$mm}{"chromosome"}}, $loc[0]-1, $loc[1]-$loc[0]+1);
      $miseq2 = substr($fas{$miR{$mm}{"chromosome"}}, $loc[2]-1, $loc[3]-$loc[2]+1) if($#loc > 2);
      if($miR{$mm}{"strand"} eq "+"){
        if($#loc > 2){
          print OUTPUT3 "$refseq\t$mm\n"."."x($loc[0]-$lim1-1).$miseq1."."x($loc[2]-$loc[1]-1).$miseq2."."x($lim2-$loc[3]+1)."\tREF\n";
        }else{
          print OUTPUT3 "$refseq\t$mm\n"."."x($loc[0]-$lim1-1).$miseq1."."x($lim2-$loc[1]+1)."\tREF\n";
        }
      }else{
        $refseq = Function->revcomp($refseq);
        $miseq1 = Function->revcomp($miseq1);
        $miseq2 = Function->revcomp($miseq2) if($#loc > 2);
        if($#loc > 2){
          print OUTPUT3 "$refseq\t$mm\n"."."x($lim2-$loc[3]+1).$miseq2."."x($loc[2]-$loc[1]-1).$miseq1."."x($loc[0]-$lim1-1)."\tREF\n";
        }else{
          print OUTPUT3 "$refseq\t$mm\n"."."x($lim2-$loc[1]+1).$miseq1."."x($loc[0]-$lim1-1)."\tREF\n";
        }
      }
      foreach my $read (sort keys %{$miR{$mm}{"read"}}){
        my @features = split /_/, $read;
        $features[2] = "" if(!$features[2]);
        if($miR{$mm}{"strand"} eq "+"){
          $seq = substr($fas{$miR{$mm}{"chromosome"}}, $features[0], $features[1]-$features[0]);
          $flank1 = $features[0] - $lim1;
          $flank2 = $lim2 - $features[1] - length($features[2]);
        }else{
          $seq = substr($fas{$miR{$mm}{"chromosome"}}, $features[0], $features[1]-$features[0]);
          $seq = Function->revcomp($seq);
          $flank1 = $lim2 - $features[1] + 1;
          $flank2 = $features[0] - $lim1 - length($features[2]) - 1;
        }
        print OUTPUT3 "."x$flank1.$seq.$features[2]."."x($flank2+1)."\t$miR{$mm}{read}{$read}{count}\n";
      }
    }
    close OUTPUT3;
    unlink "$tag.out2";
  }
  my $par = join " ", @par;
  system ("Rscript --vanilla ".$prefix."/scripts/bubble_plot.R ".$par);
}

1;
__END__
