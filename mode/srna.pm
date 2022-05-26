#!/usr/bin/env perl
package mode::srna;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Path qw/remove_tree/;
use Function;
use Ref;
use validate_options;
use input;
use Cwd 'abs_path';

command_short_description q[Analysis for small RNA-seq];
command_long_description q[Analysis for small RNA-seq];
command_usage q[pRNASeqTools srna [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT]=[file1]+[file2] ... ];

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
  default => 'rRNA,total',
  documentation => q[Method for normalization, seperated by comma. Default: rRNA, total.],
);
option 'binsize' => (
  is => 'rw',
  isa => 'Num',
  default => 100,
  documentation => q[Windows size to slice the genome.],
);
option 'mode' => (
  is => 'rw',
  isa => 'Str',
  default => 'bulk',
  documentation => q[method for assigning multiple mapped reads. Allowed: bulk, sc],
);
option 'pattern' => (
  is => 'rw',
  isa => 'Str',
  default => 'NNNNNNNNCA',
  documentation => q[UMI pattern, only effective with mode 'sc'],
);
option 'promoter' => (
  is => 'rw',
  isa => 'Num',
  default => 1000,
  documentation => q[Length of gene promoter],
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
  my $norm = $options{'norm'};
  my @norms = split /,/, $norm;
  my $mask = $options{'mask'};
  my $mode = $options{'mode'};
  my $pattern = $options{'pattern'};
  my $promoterLength = $options{'promoter'};

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
    Ref->splitgff($prefix, $genome, $promoterLength);
    if(defined $mask){
      if($mask =~ /^~\/(.+)/){
        $mask = $ENV{"HOME"}."/".$1;
      }elsif($mask !~ /^\//){
        $mask = abs_path "../".$mask;
      }
      symlink $mask, "mask.fa";
      system ("bowtie-build -q mask.fa mask");
    }
    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];

      print $main::tee "\nMapping $tag...\n";

      $file = Function->SRR($file, $thread);
      Function->unzip($file, $tag);
      if($mode eq 'sc'){

        print $main::tee "\nTrimming $tag...\n";

        system ("umi_tools extract -p ".$pattern." -I ".$tag.".fastq -S ".$tag.".fq");
        if(defined $adaptor){
          system ("cutadapt -j ".$thread." -m 18 -M 42 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fq 2>&1");
        }
        open DUP, "<".$tag."_trimmed.fastq" or die $!;
        my (%dedup, $outSeq, @seqNames, @outSeqNames);
        while (my $dup = <DUP>){
          chomp $dup;
          if($. % 4 == 0){
            $dedup{$outSeqNames[2]}{$outSeq}{"quality"} = $dup;
          }elsif($. % 4 == 2){
            $outSeq = $dup;
            $dedup{$outSeqNames[2]}{$outSeq}{"id"} = $outSeqNames[0];
            $dedup{$outSeqNames[2]}{$outSeq}{"count"} ++;
          }elsif($. % 4 == 1){
            @seqNames = split / /, $dup;
            @outSeqNames = split /_/, $seqNames[0];
          }
        }
        close DUP;
        open DEDUP, ">".$tag.".fastq" or die $!;
        foreach my $umi (keys %dedup){
          foreach my $outputSeq (keys %{$dedup{$umi}}){
            print DEDUP "$dedup{$umi}{$outputSeq}{id}\_$umi\_$dedup{$umi}{$outputSeq}{count} $seqNames[1]\n$outputSeq\n+\n$dedup{$umi}{$outputSeq}{quality}\n";
          }
        }
        close DEDUP;
        unlink $tag."_trimmed.fastq", $tag.".fq";
      }elsif($mode eq 'bulk'){
        if(defined $adaptor){

          print $main::tee "\nTrimming $tag...\n";

          system ("cutadapt -j ".$thread." -m 18 -M 42 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
          rename $tag."_trimmed.fastq", $tag.".fastq";
        }
      }
      if(defined $mask){
        system ("bowtie -v 0 -a --un tmp.fastq -p ".$thread." -t mask ".$tag.".fastq ".$tag.".mask.out 2>&1");
        rename "tmp.fastq", $tag.".fastq";
        unlink $tag.".mask.out";
      }

      print $main::tee "\nStart mapping...\n";

      system ("bowtie -v 2 -a -p ".$thread." -t ".$prefix."/reference/lsu_rrna ".$tag.".fastq ".$tag.".rRNA.out 2>&1");
      system ("awk -F \"\t\" \'BEGIN{x=0;y=0;z=0}{if(\$2==\"+\"){if(/$genome\_LSU/){x++};if(/$genome\_SSU/){y++};if(/$genome\_U6/){z++}}}END{print \"rRNA\t\"x\"\\nSSU\t\"y\"\\nU6\t\"z}\' ".$tag.".rRNA.out > ".$tag.".nf");

      system ("ShortStack --outdir ShortStack_".$tag." --align_only --bowtie_m 1000 --ranmax 50 --mmap ".$mmap." --mismatches 0 --bowtie_cores ".$thread." --nohp --readfile ".$tag.".fastq --genomefile ".$prefix."/reference/".$genome."_chr_all.fasta 2>&1");

      print $main::tee "\nAlignment Completed!\n";

      system ("samtools view -h ShortStack_".$tag."/".$tag.".bam > ".$tag);
      system ("awk '{if(\$0~/^@/) print > (FILENAME\".unmapped.sam\"); if(\$10!=\"*\" && \$3!=\"*\") print > (FILENAME\".sam\"); if(\$10!=\"*\" && \$3==\"*\") print > (FILENAME\".unmapped.sam\")}' ".$tag);
      system ("samtools view -Sb ".$tag.".unmapped.sam > ".$tag.".unmapped.bam");
      system ("samtools view -Sb ".$tag.".sam > ".$tag.".bam");

      unlink ($tag, $tag.".unmapped.sam", $tag.".sam", $tag.".fastq", $tag.".rRNA.out");
      remove_tree "ShortStack_".$tag;

      print $main::tee "\nConverting BAM to BED...\n";

      system ("bamToBed -bed12 -i ".$tag.".bam > ".$tag.".bed");

      print $main::tee "\nGenerating individual files...\n";

      system ("awk -F \"\t\" '{a=substr(FILENAME,1,length(FILENAME)-3);if(\$11>=18 && \$11 <= 26) {print \$0 > (a\$11\".bed\")}}' ".$tag.".bed");
      system ("awk '!a[\$4]++' ".$tag.".bed \| awk '{print \$11}' \| sort \| uniq -c \| awk '{OFS=\"\t\"; print \$2, \$1}' > ".$tag.".len_dist.txt");
      system ("awk '{if(\$1<=28){n+=\$2}}END{print \"total\t\"n}' ".$tag.".len_dist.txt >> ".$tag.".nf");

      print $main::tee "\nLength distribution summary done!\nCounting start...\n";

      open NF, $tag.".nf" or die $!;
      my %normhash;
      while (my $line = <NF>){
        chomp $line;
        my ($mnorm, $rc) = split /\t/, $line;
        $normhash{$mnorm} = $rc;
      }
      close NF;
      foreach my $mnorm (@norms){
        if(exists $normhash{$mnorm}){
          mode::srna->count($mnorm, $normhash{$mnorm}, $prefix, $genome, $binsize, $tag);
        }
      }

      print $main::tee "Counting Completed!\n";

      unlink glob ($tag."*.bed");
    }
    unlink glob ("mask*") if(defined $mask);
    if(!$mappingonly && $#par > 1){
      foreach my $mnorm (@norms){
        mode::srna->sta($mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $promoterLength, $par);
      }
    }
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
    foreach my $mnorm (@norms){
      mode::srna->sta($mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $promoterLength, $par);
    }
    unlink glob ("*.count"), glob ("*.nf");
  }
  unlink glob ("*.gff");
}


sub sta {
  my ($self, $mnorm, $prefix, $genome, $foldchange, $pvalue, $binsize, $promoterLength, $par) = @_;

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
    system ("bedGraphToBigWig ".$bg." ".$prefix."/reference/".$genome."_chr_all.fasta.fai ".$bw);
    unlink $bg;
  }

  my %ann = Ref->ann($prefix, $genome, $binsize, $promoterLength);
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

	print $main::tee "\nDE miRNA analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DEM.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
	print $main::tee "\nDS gene analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DSG.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
	print $main::tee "\nDS TE analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DST.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
	print $main::tee "\nDS Promoter analysis...\nNormalization $mnorm\tFold Change $foldchange\tP Value $pvalue\n";
	system ("Rscript --vanilla ".$prefix."/scripts/DSP.R ".$mnorm." ".$pvalue." ".$foldchange." ".$par);
}

sub count {
  my ($self, $mnorm, $rc, $prefix, $genome, $binsize, $tag) = @_;
  my (%mirna, $miseq, %count);

  system ("bedtools intersect -a ".$prefix."/reference/".$genome."_miRNA_miRNA_star.gff -b ".$tag.".bam -wa -f 1 -r -c \|awk '{print \$9\"\t\"\$10}' > ".$tag.".miRNA.tmp");
  my %fas = Ref->fas($prefix, $genome);
  open MIA, "<$tag".".miRNA.tmp" or die $!;
  while (my $mi = <MIA>){
    chomp $mi;
    my @row = split /\t/, $mi;
    $count{$row[0]} = $row[1];
  }
  close MIA;
  unlink $tag.".miRNA.tmp";
  open MI, "$prefix/reference/".$genome."_miRNA_miRNA_star.gff" or die $!;
  while (my $mi = <MI>){
    chomp $mi;
    my @row = split /\t/, $mi;
    $miseq = substr($fas{$row[0]}, $row[3]-1, $row[4]-$row[3]+1);
    $miseq = Function->revcomp($miseq) if($row[6] eq "-");
    $mirna{$miseq}{"name"} .= $row[8].";";
    $mirna{$miseq}{"count"} += $count{$row[8]};
  }
  close MI;
  open MIO, ">".$tag.".miRNA.annotated.count" or die $!;
  foreach $miseq (sort keys %mirna){
    chop $mirna{$miseq}{"name"};
    print MIO "$mirna{$miseq}{name}\t$mirna{$miseq}{count}\n";
  }
  close MIO;
  %count = ();
  %mirna = ();
  opendir my $dir, "." or die $!;
  my @dir = grep {/$tag.*[0-9]\.bed/} readdir $dir;
  close $dir;
  foreach my $sbed (@dir){
    my $nrc = 1000000 / $rc;

    (my $bgp = $sbed) =~ s/bed$/p.bedgraph/;
    system ("bedtools genomecov -split -strand + -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bgp);
    (my $bgpo = $bgp) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwpo = $bgp) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -strand + -scale ".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bgpo);
    system ("bedGraphToBigWig ".$bgpo." ".$prefix."/reference/".$genome."_chr_all.fasta.fai ".$bwpo);

    (my $bgn = $sbed) =~ s/bed$/n.bedgraph/;
    system ("bedtools genomecov -split -strand - -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bgn);
    (my $bgno = $bgn) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwno = $bgn) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -strand - -scale -".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bgno);
    system ("bedGraphToBigWig ".$bgno." ".$prefix."/reference/".$genome."_chr_all.fasta.fai ".$bwno);

    (my $bg = $sbed) =~ s/bed$/bedgraph/;
    system ("bedtools genomecov -split -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bg);
    (my $bgo = $bg) =~ s/bedgraph$/$mnorm\.bedgraph/;
    (my $bwo = $bg) =~ s/bedgraph$/$mnorm\.bw/;
    system ("bedtools genomecov -split -scale ".$nrc." -bg -i ".$sbed." -g ".$prefix."/reference/".$genome."_chr_all.fasta.fai > ".$bgo);
    system ("bedGraphToBigWig ".$bgo." ".$prefix."/reference/".$genome."_chr_all.fasta.fai ".$bwo);

    unlink($bgp, $bgn, $bg, $bgpo, $bgno, $bgo);

    #Counting miRNA

    system ("bedtools intersect -a ".$prefix."/reference/".$genome."_miRNA_miRNA_star.gff -b ".$sbed." -wa -f 0.95 -c \|awk -v x=".$rc." '{print \$9\"\t\"\$10 * 1000000 / x}' > ".$sbed.".tmp");
    open TMP, "<$sbed.tmp" or die $!;
    while (my $line = <TMP>){
      chomp $line;
      my @row = split /\t/, $line;
      $count{"mir"}{$row[0]}{$sbed}{"n"} = $row[1];
    }
    close TMP;
    system ("bedtools intersect -a ".$prefix."/reference/".$genome."_miRNA_miRNA_star.gff -b ".$sbed." -wa -f 0.95 -c \|awk '{print \$9\"\t\"\$10}' > ".$sbed.".tmp");
    open TMP, "<$sbed.tmp" or die $!;
    while (my $line = <TMP>){
      chomp $line;
      my @row = split /\t/, $line;
      $count{"mir"}{$row[0]}{$sbed}{"r"} = $row[1];
    }
    close TMP;
    unlink $sbed.".tmp";

    #Counting small RNAs in binsize regions

    open BED, "<$sbed" or die $!;
    while (my $bed = <BED>){
      chomp $bed;
      my @row = split /\t/, $bed;
      my $bin = int(($row[1]+1)/$binsize);
      $count{"r100"}{$row[0]}{$bin}{$sbed}{"r"} ++;
    }
    close BED;

    my %length = Ref->lengthofchrom($prefix, $genome, $binsize);
    foreach my $chr (sort keys %length){
        for(my $bi=0;$bi<=$length{$chr};$bi++){
        if(exists $count{"r100"}{$chr}{$bi}{$sbed}{"r"}){
          $count{"r100"}{$chr}{$bi}{$sbed}{"n"} = $count{"r100"}{$chr}{$bi}{$sbed}{"r"} * 1000000 / $rc;
        }else{
          $count{"r100"}{$chr}{$bi}{$sbed}{"n"} = 0;
          $count{"r100"}{$chr}{$bi}{$sbed}{"r"} = 0;
        }
      }
    }

    #Counting small RNAs in genome features

    system ("bedtools intersect -a gene.gff -b ".$sbed." -wa -c > ".$sbed.".gene.tmp");
    open FILE, "$sbed.gene.tmp" or die $!;
    while(my $cc = <FILE>){
      chomp $cc;
      my @row = split /\t/, $cc;
      $count{"gene"}{$row[8]}{$sbed}{"r"} = $row[9];
      $count{"gene"}{$row[8]}{$sbed}{"n"} = $row[9] * 1000000 / $rc;
    }
    close FILE;
    unlink $sbed.".gene.tmp";

    system ("bedtools intersect -a te.gff -b ".$sbed." -wa -c > ".$sbed.".te.tmp");
    open FILE, "$sbed.te.tmp" or die $!;
    while(my $dd = <FILE>){
      chomp $dd;
      my @row = split /\t/, $dd;
      $count{"te"}{$row[8]}{$sbed}{"r"} = $row[9];
      $count{"te"}{$row[8]}{$sbed}{"n"} = $row[9] * 1000000 / $rc;
    }
    close FILE;
    unlink $sbed.".te.tmp";

    system ("bedtools intersect -a promoter.gff -b ".$sbed." -wa -c > ".$sbed.".promoter.tmp");
    open FILE, "$sbed.promoter.tmp" or die $!;
    while(my $ee = <FILE>){
      chomp $ee;
      my @row = split /\t/, $ee;
      $count{"promoter"}{$row[8]}{$sbed}{"r"} = $row[9];
      $count{"promoter"}{$row[8]}{$sbed}{"n"} = $row[9] * 1000000 / $rc;
    }
    close FILE;
    unlink $sbed.".promoter.tmp";
  }

  #Printing counting results

  open GEN1, ">$tag.count" or die $!;
  open GEN2, ">$tag.$mnorm.norm.count" or die $!;
  foreach my $chr (sort keys %{$count{"r100"}}){
    foreach my $bi (sort {$a <=> $b} keys %{$count{"r100"}{$chr}}){
      print GEN1 "$chr\_$bi";
      print GEN2 "$chr\_$bi";
      foreach my $leng (sort keys %{$count{"r100"}{$chr}{$bi}}){
        print GEN1 "\t$count{r100}{$chr}{$bi}{$leng}{r}";
        print GEN2 "\t$count{r100}{$chr}{$bi}{$leng}{n}";
      }
      print GEN1 "\n";
      print GEN2 "\n";
    }
  }
  close GEN1;
  close GEN2;

  open GEN1, ">$tag.gene.count" or die $!;
  open GEN2, ">$tag.gene.$mnorm.norm.count" or die $!;
  foreach my $gene (sort keys %{$count{"gene"}}){
    print GEN1 "$gene";
    print GEN2 "$gene";
    foreach my $leng (sort keys %{$count{"gene"}{$gene}}){
      print GEN1 "\t$count{gene}{$gene}{$leng}{r}";
      print GEN2 "\t$count{gene}{$gene}{$leng}{n}";
    }
    print GEN1 "\n";
    print GEN2 "\n";
  }
  close GEN1;
  close GEN2;

  open GEN1, ">$tag.TE.count" or die $!;
  open GEN2, ">$tag.TE.$mnorm.norm.count" or die $!;
  foreach my $gene (sort keys %{$count{"te"}}){
    print GEN1 "$gene";
    print GEN2 "$gene";
    foreach my $leng (sort keys %{$count{"te"}{$gene}}){
      print GEN1 "\t$count{te}{$gene}{$leng}{r}";
      print GEN2 "\t$count{te}{$gene}{$leng}{n}";
    }
    print GEN1 "\n";
    print GEN2 "\n";
  }
  close GEN1;
  close GEN2;

  open GEN1, ">$tag.promoter.count" or die $!;
  open GEN2, ">$tag.promoter.$mnorm.norm.count" or die $!;
  foreach my $gene (sort keys %{$count{"promoter"}}){
    print GEN1 "$gene";
    print GEN2 "$gene";
    foreach my $leng (sort keys %{$count{"promoter"}{$gene}}){
      print GEN1 "\t$count{promoter}{$gene}{$leng}{r}";
      print GEN2 "\t$count{promoter}{$gene}{$leng}{n}";
    }
    print GEN1 "\n";
    print GEN2 "\n";
  }
  close GEN1;
  close GEN2;

  #miRNA

  open MI, "$prefix/reference/".$genome."_miRNA_miRNA_star.gff" or die $!;
  while (my $mi = <MI>){
    chomp $mi;
    my @row = split /\t/, $mi;
    $miseq = substr($fas{$row[0]}, $row[3]-1, $row[4]-$row[3]+1);
    $miseq = Function->revcomp($miseq) if($row[6] eq "-");
    $mirna{$miseq}{"name"} .= $row[8].";";
    foreach my $leng (keys %{$count{"mir"}{$row[8]}}){
      $mirna{$miseq}{$leng}{"r"} += $count{"mir"}{$row[8]}{$leng}{"r"};
      $mirna{$miseq}{$leng}{"n"} += $count{"mir"}{$row[8]}{$leng}{"n"};
    }
  }
  close MI;

  open MI1, ">$tag.miRNA.count" or die $!;
  open MI2, ">$tag.miRNA.$mnorm.norm.count" or die $!;
  foreach $miseq (sort keys %mirna){
    chop $mirna{$miseq}{"name"};
    print MI1 "$mirna{$miseq}{name}";
    print MI2 "$mirna{$miseq}{name}";
    foreach my $leng (sort keys %{$mirna{$miseq}}){
      next if($leng eq "name");
      print MI1 "\t$mirna{$miseq}{$leng}{r}";
      print MI2 "\t$mirna{$miseq}{$leng}{n}";
    }
    print MI1 "\n";
    print MI2 "\n";
  }
  close MI1;
  close MI2;
}

1;

__END__
