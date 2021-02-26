#!/usr/bin/env perl
package mode::risi;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Path qw/remove_tree/;
use Function;
use Ref;
use validate_options;
use input;
use Cwd 'abs_path';

command_short_description q[Analysis for risiRNA];
command_long_description q[Analysis for risiRNA];
command_usage q[pRNASeqTools risi [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT]=[file1]+[file2] ... ];

option 'mmap' => (
  is => 'rw',
  isa => 'Str',
  default => 'u',
  documentation => q[method for assigning multiple mapped reads. Allowed: u, n, f, r],
);
option 'no-mapping' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Just perform the statistic analysis],
);
option 'foldchange' => (
  is => 'rw',
  isa => 'Num',
  default => '1.5',
  documentation => q[Threshold for DEG in fold change.],
);
option 'pvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '0.01',
  documentation => q[Threshold for DEG in P value. Default: 0.01],
);
option 'mapping-only' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Do not perform the statistic analysis],
);
option 'norm' => (
  is => 'rw',
  isa => 'Str',
  default => 'total',
  documentation => q[Method for normalization, seperated by comma. Default: rRNA, total.],
);
option 'binsize' => (
  is => 'rw',
  isa => 'Num',
  default => 10,
  documentation => q[Windows size to slice the genome.],
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
  my $nomapping = $options{'no-mapping'};
  my $mappingonly = $options{'mapping-only'};
  my $foldchange = $options{'foldchange'};
  my $pvalue = $options{'pvalue'};
  my $binsize = $options{'binsize'};
  my $pattern = $options{'pattern'};
  my $mnorm;
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
    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];
      print $main::tee "\nMapping $tag...\n";
      $file = Function->SRR($file, $thread);
      Function->unzip($file, $tag);
      if(defined $adaptor){
        print $main::tee "\nTrimming $tag...\n";
        system ("cutadapt -j ".$thread." -m 18 -M 42 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
        rename $tag."_trimmed.fastq", $tag.".fastq";
      }
      print $main::tee "\nStart mapping...\n";
      system ("ShortStack --outdir ShortStack_".$tag." --align_only --bowtie_m 1000 --ranmax 50 --mmap ".$mmap." --mismatches 1 --bowtie_cores ".$thread." --nohp --readfile ".$tag.".fastq --genomefile ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta 2>&1");
      print $main::tee "\nAlignment Completed!\n";
      system ("samtools view -h ShortStack_".$tag."/".$tag.".bam > ".$tag);
      system ("awk '{if(\$0~/^@/) print > (FILENAME\".unmapped.sam\"); if(\$10!=\"*\" && \$3!=\"*\") print > (FILENAME\".sam\"); if(\$10!=\"*\" && \$3==\"*\") print > (FILENAME\".unmapped.sam\")}' ".$tag);
      system ("samtools view -Sb ".$tag.".unmapped.sam > ".$tag.".unmapped.bam");
      system ("samtools view -Sb ".$tag.".sam > ".$tag.".bam");
      unlink ($tag, $tag.".unmapped.sam", $tag.".sam", $tag.".fastq");
      remove_tree "ShortStack_".$tag;
      print $main::tee "\nConverting BAM to BED...\n";
      system ("bamToBed -bed12 -i ".$tag.".bam > ".$tag.".bed");
      print $main::tee "\nGenerating individual files...\n";
      system ("awk -F \"\t\" '{a=substr(FILENAME,1,length(FILENAME)-3);if(\$11>=18 && \$11 <= 26) {print \$0 > (a\$11\".bed\")}}' ".$tag.".bed");
      system ("awk '!a[\$4]++' ".$tag.".bed \| awk '{print \$11}' \| sort \| uniq -c \| awk '{OFS=\"\t\"; print \$2, \$1}' > ".$tag.".len_dist.txt");
      system ("awk '{n+=\$2}END{print \"total\t\"n}' ".$tag.".len_dist.txt >> ".$tag.".nf");
      print $main::tee "\nLength distribution summary done!\nCounting start...\n";
      open NF, $tag.".nf" or die $!;
      my (%normhash, $rc);
      while (my $line = <NF>){
        chomp $line;
        ($mnorm, $rc) = split /\t/, $line;
        $normhash{$mnorm} = $rc;
      }
      close NF;
      mode::risi->count($mnorm, $normhash{$mnorm}, $prefix, $genome, $binsize, $tag) if(exists $normhash{$mnorm});
      print $main::tee "Counting Completed!\n";
      unlink glob ($tag."*.bed");
    }
    mode::risi->sta($mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $par) if(!$mappingonly && $#par > 1);
  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".nf", $pre.".nf" or die $!;
      opendir my $originalFolder, "../" or die $!;
      my @fileList = readdir $originalFolder;
      foreach my $file (@fileList){
        symlink "../".$file, $file if($file =~ /^$pre.+count$/ and $file !~ /norm/);
      }
      closedir $originalFolder;
    }
    mode::risi->sta($mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $par);
    unlink glob ("*.count"), glob ("*.nf");
  }
  unlink glob ("*.gff");
}


sub sta {
  my ($self, $mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $par) = @_;
	print $main::tee "\nDSR analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DSR.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
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
    (my $bw = $hcsv) =~ s/csv$/bw/;
    system ("bedGraphToBigWig ".$bg." ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai ".$bw);
    unlink $bg;
  }
	print $main::tee "\nDS feature analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DSF.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
}
sub count {
  my ($self, $mnorm, $rc, $prefix, $genome, $binsize, $tag) = @_;
  my %count;
  my %length = Ref->lengthofchrom($prefix, $genome."_rDNA", $binsize);
  opendir my $dir, "." or die $!;
  my @dir = grep {/$tag.*[0-9]\.bed/} readdir $dir;
  close $dir;
  foreach my $sbed (@dir){
    my $nrc = 1000000 / $rc;
    (my $bgp = $sbed) =~ s/bed$/p.bedgraph/;
    system ("bedtools genomecov -split -strand + -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bgp);
    (my $bgpo = $bgp) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwpo = $bgp) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -strand + -scale ".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bgpo);
    system ("bedGraphToBigWig ".$bgpo." ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai ".$bwpo);
    (my $bgn = $sbed) =~ s/bed$/n.bedgraph/;
    system ("bedtools genomecov -split -strand - -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bgn);
    (my $bgno = $bgn) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwno = $bgn) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -strand - -scale -".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bgno);
    system ("bedGraphToBigWig ".$bgno." ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai ".$bwno);
    (my $bg = $sbed) =~ s/bed$/bedgraph/;
    system ("bedtools genomecov -split -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bg);
    (my $bgo = $bg) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwo = $bg) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -scale ".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai > ".$bgo);
    system ("bedGraphToBigWig ".$bgo." ".$prefix."/reference/".$genome."_rDNA_chr_all.fasta.fai ".$bwo);
    unlink($bgp, $bgn, $bg, $bgpo, $bgno, $bgo);
    #Counting small RNAs in binsize regions
    open BED, "<$sbed" or die $!;
    my $leng = $1 if($sbed =~ /(\d+)\.bed/);
    while (my $bed = <BED>){
      chomp $bed;
      my @row = split /\t/, $bed;
      my $bin = int(($row[1]+1)/$binsize);
      $count{"r10"}{$row[0]}{$bin}{$row[5]}{$leng} ++;
    }
    close BED;
    #Counting small RNAs in rDNA features
    system ("bedtools intersect -a ".$sbed." -b ".$prefix."/reference/".$genome."_rDNA.gff -wb > ".$sbed.".tmp");
    my %tmpHash;
    open FILE, "$sbed.tmp" or die $!;
    while(my $cc = <FILE>){
      chomp $cc;
      my @row = split /\t/, $cc;
      $tmpHash{$row[3]}{"strand"} = $row[5];
      $tmpHash{$row[3]}{"feature"} .= $row[20].",";
      $tmpHash{$row[3]}{"length"} = $row[10];
    }
    close FILE;
    foreach my $read (sort keys %tmpHash){
      chop $tmpHash{$read}{"feature"};
      $count{"feature"}{$tmpHash{$read}{feature}}{$tmpHash{$read}{strand}}{$tmpHash{$read}{length}} ++;
    }
    unlink $sbed.".tmp";
  }
  #Printing counting results
  my $tmpNorm;
  open GEN1, ">$tag.count" or die $!;
  open GEN2, ">$tag.$mnorm.norm.count" or die $!;
  foreach my $chr (sort keys %{$count{"r10"}}){
    for (my $bi=0;$bi<=$length{$chr};$bi++){
      foreach my $strand (("+","-")){
        print GEN1 "$chr\_$bi\_$strand";
        print GEN2 "$chr\_$bi\_$strand";
        for(my $leng=18;$leng<=26;$leng++){
          if(exists $count{"r10"}{$chr}{$bi}{$strand}{$leng}){
            print GEN1 "\t$count{r10}{$chr}{$bi}{$strand}{$leng}";
            $tmpNorm = $count{"r10"}{$chr}{$bi}{$strand}{$leng} * 1000000 / $rc;
            print GEN2 "\t$tmpNorm";
          }else{
            print GEN1 "\t0";
            print GEN2 "\t0";
          }
        }
        print GEN1 "\n";
        print GEN2 "\n";
      }
    }
  }
  close GEN1;
  close GEN2;
  open GEN1, ">$tag.feature.count" or die $!;
  open GEN2, ">$tag.feature.$mnorm.norm.count" or die $!;
  foreach my $feature (sort keys %{$count{"feature"}}){
    foreach my $strand (("+","-")){
      print GEN1 "$feature\_$strand";
      print GEN2 "$feature\_$strand";
      for(my $leng=18;$leng<=26;$leng++){
        if(exists $count{"feature"}{$feature}{$strand}{$leng}){
          print GEN1 "\t$count{feature}{$feature}{$strand}{$leng}";
          $tmpNorm = $count{"feature"}{$feature}{$strand}{$leng} * 1000000 / $rc;
          print GEN2 "\t$tmpNorm";
        }else{
          print GEN1 "\t0";
          print GEN2 "\t0";
        }
      }
      print GEN1 "\n";
      print GEN2 "\n";
    }
  }
  close GEN1;
  close GEN2;
}

1;

__END__
