#!/usr/bin/env perl
package mode::mrna;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Path qw/remove_tree make_path/;
use File::Cat;
use Function;
use Ref;
use validate_options;
use input;
use Cwd 'abs_path';

command_short_description q[Analysis for RNA-seq];
command_long_description q[Analysis for RNA-seq];
command_usage q[pRNASeqTools mrna [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT1]=[file1]+[file2] ... --treatment [TREATMENT2]=[file1]+[file2] ... ];

option 'total' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Perform the total RNA-seq analysis (including ncRNA)],
);
option 'mode' => (
  is => 'rw',
  isa => 'Num',
  default => '1',
  documentation => q[Specify the analyzing mode, allowed: 1 for input file in in fastq format, 2 for mapping and count only, 3 for input file is in bam format, 4 for input file is tab-delimited readcount table],
);
option 'foldchange' => (
  is => 'rw',
  isa => 'Num',
  default => '2',
  documentation => q[Threshold for DEG in fold change.],
);
option 'pvalue' => (
  is => 'rw',
  isa => 'Num',
  default => '0.01',
  documentation => q[Threshold for DEG in P value.],
);
option 'fdr' => (
  is => 'rw',
  isa => 'Num',
  default => '1',
  documentation => q[Threshold for DEG in FDR. Will overwrite P value threshold],
);
option 'DESeq2Norm' => (
  is => 'rw',
  isa => 'Str',
  default => 'DESeq2',
  documentation => q[Normalization method for DESeq2, "DESeq2" and "RPM" are allowed],
);
option 'seqStrategy' => (
  is => 'rw',
  isa => 'Str',
  documentation => q[Sequencing stragety, must specify performing bam mode. Allowed value: paired, single]
);
option 'genomeSize' => (
  is => 'rw',
  isa => 'Num',
  default => '10',
  documentation => q[Parameter 'genomeSAindexNbases' for STAR. Default: 10],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $foldchange = $options{'foldchange'};
  my $pvalue = $options{'pvalue'};
  my $fdr = $options{'fdr'};
  my $mode = $options{'mode'};
  my $control = $options{'control'};
  my $treatment = $options{'treatment'};
  my $norm = $options{'DESeq2Norm'};
  my $mask = $options{'mask'};
  my $total = $options{'total'};
  my $seqStrategy = $options{'seqStrategy'};
  my $genomeSize = $options{'genomeSize'};

  my $mapping = 1,;
  my $count = 1;
  my $de = 1;
  $de = 0 if($mode == 2);
  $mapping = 0 if($mode == 3);
  if($mode == 4){
    $mapping = 0;
    $count = 0;
  }
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

  if($mapping){
    if(defined $mask){
      if($mask =~ /^~\/(.+)/){
        $mask = $ENV{"HOME"}."/".$1;
      }elsif($mask !~ /^\//){
        $mask = abs_path "../".$mask;
      }
      symlink $mask, "mask.fa";
      system ("bowtie-build -q mask.fa mask");
    }

    print $main::tee "\nBuilding STAR genome index ...\n";

    remove_tree "Genome" if(-e "Genome");
    make_path "Genome";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeSAindexNbases ".$genomeSize." --genomeFastaFiles ".$prefix."/reference/".$genome."_chr_all.fasta --sjdbGTFfile ".$prefix."/reference/".$genome."_genes.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000");
    if($total){
      system ("gffread -T -o ".$genome."_genes.gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    }else{
      system ("gffread -T -C -o ".$genome."_genes.gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    }
    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];
  		print $main::tee "\nProcessing $tag...\n";
  		if($file !~ /,/){
        $seqStrategy = "single";
			  my @files = Function->SRR($file, $thread);
        if($#files == 0){
      	  Function->unzip($files[0], $tag);
      	  if(defined $adaptor){
				 	  system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
          	rename $tag."_trimmed.fastq", $tag.".fastq";
        	}
          if(defined $mask){
            system ("bowtie -v 0 -a --un tmp.fastq -p ".$thread." -t mask ".$tag.".fastq ".$tag.".mask.out 2>&1");
            rename "tmp.fastq", $tag.".fastq";
            unlink $tag.".mask.out";
          }
    			system ("STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
    			unlink ($tag.".fastq");
        }else{
          $seqStrategy = "paired";
          Function->unzip($files[0], $tag."_R1");
          Function->unzip($files[1], $tag."_R2");
          if(defined $adaptor){
            system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
            rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
            rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
          }
          if(defined $mask){
            system ("bowtie -v 0 -a --un tmp.fastq -p ".$thread." -t mask -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq ".$tag.".mask.out 2>&1");
            rename "tmp_1.fastq", $tag."_R1.fastq";
            rename "tmp_2.fastq", $tag."_R2.fastq";
            unlink $tag.".mask.out";
          }
          system ("STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          unlink ($tag."_R1.fastq", $tag."_R2.fastq");
        }
  		}else{
       	my ($file1, $file2) = split /,/, $file;
        $seqStrategy = "paired";
				Function->unzip($file1, $tag."_R1");
				Function->unzip($file2, $tag."_R2");
  		  if(defined $adaptor){
          system ("cutadapt -j ".$thread." -m 20 --trim-n -a ".$adaptor." -A ".$adaptor." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
          rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
        }
        if(defined $mask){
          system ("bowtie -v 0 -a --un tmp.fastq -p ".$thread." -t mask -1 ".$tag."_R1.fastq -2 ".$tag."_R2.fastq ".$tag.".mask.out 2>&1");
          rename "tmp_1.fastq", $tag."_R1.fastq";
          rename "tmp_2.fastq", $tag."_R2.fastq";
          unlink $tag.".mask.out";
        }
        print $main::tee "\nMapping...\n";
        system ("STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
  			unlink ($tag."_R1.fastq", $tag."_R2.fastq");
  		}
  		rename "Aligned.sortedByCoord.out.bam", $tag.".bam";
  		cat 'Log.final.out', \*STDOUT;
      mode::mrna->countBAM($tag, $thread, $seqStrategy, $prefix, $genome);
    }
    unlink ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab", $genome."_genes.gtf");
    unlink glob ("total.count*");
    unlink glob ("mask*") if(defined $mask);
    remove_tree "Genome";
    if($de){
      print $main::tee "\nFinding DEG...\nFold Change\t$foldchange\tP Value\t$pvalue\tFDR\t$fdr\n";
      system ("Rscript --vanilla ".$prefix."/scripts/DEG.R ".$norm." ".$pvalue." ".$fdr." ".$foldchange." ".$prefix." ".$genome." ".$par);
      opendir my $dir, "." or die $!;
      my @dir = grep {/csv$/} readdir $dir;
      closedir $dir;
      my %gann = Ref->gann($prefix, $genome);
      foreach my $csv (@dir){
        open CSV, "$csv" or die $!;
    		open TMP, ">tmp" or die $!;
    		while(my $bb = <CSV>){
    			chomp $bb;
    			my @row = split /,/, $bb;
    			$row[0] =~ s/"//g;
    			if(exists $gann{$row[0]}){
    				print TMP "$bb,$gann{$row[0]}\n";
    			}else{
    				print TMP "$bb,Mapman,type,short,description,long\n";
    			}
    		}
    		close CSV;
    		close TMP;
    		rename "tmp", $csv;
      }
	  }
  }else{
    foreach my $pre (@tags){
      if($count){
        symlink "../".$pre.".bam", $pre.".bam" or die $!;
        if($total){
          system ("gffread -T -o ".$genome."_genes.gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
        }else{
          system ("gffread -T -C -o ".$genome."_genes.gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
        }
        mode::mrna->countBAM($pre, $thread, $seqStrategy, $prefix, $genome);
        unlink $genome."_genes.gtf";
      }else{
        symlink "../".$pre.".txt", $pre.".txt" or die $!;
      }
    }
    print $main::tee "\nFinding DEG...\nFold Change\t$foldchange\tP Value\t$pvalue\tFDR\t$fdr\n";
    system ("Rscript --vanilla ".$prefix."/scripts/DEG.R ".$norm." ".$pvalue." ".$fdr." ".$foldchange." ".$prefix." ".$genome." ".$par);
    opendir my $dir, "." or die $!;
    my @dir = grep {/csv$/} readdir $dir;
    closedir $dir;
    my %gann = Ref->gann($prefix, $genome);
    foreach my $csv (@dir){
      open CSV, "$csv" or die $!;
      open TMP, ">tmp" or die $!;
      while(my $bb = <CSV>){
        chomp $bb;
        my @row = split /,/, $bb;
        $row[0] =~ s/"//g;
        if(exists $gann{$row[0]}){
          print TMP "$bb,$gann{$row[0]}\n";
        }else{
          print TMP "$bb,Mapman,type,short,description,long\n";
        }
      }
      close CSV;
      close TMP;
      rename "tmp", $csv;
    }
    if($count){
      unlink glob "*_?.bam";
    }else{
      unlink glob "*_?.txt";
    }
  }
}

sub countBAM {
	my ($self, $tag, $thread, $seqStrategy, $prefix, $genome) = @_;

	system ("samtools index ".$tag.".bam");
	system ("bamCoverage -b ".$tag.".bam --skipNAs -bs 5 -p ".$thread." --minMappingQuality 10 --normalizeUsing CPM -o ".$tag.".bw");

	print $main::tee "\nStart counting...\n";

	my %count = ();
	my $countSum = 0;
	if($seqStrategy eq "single"){
		system ("featureCounts -T ".$thread." -O -G ".$prefix."/reference/".$genome."_chr_all.fasta -s 0 -a ".$genome."_genes.gtf -o total.count ".$tag.".bam 2>&1");
	}elsif($seqStrategy eq "paired"){
		system ("featureCounts -T ".$thread." -p --countReadPairs -BCO -G ".$prefix."/reference/".$genome."_chr_all.fasta -s 0 -a ".$genome."_genes.gtf -o total.count ".$tag.".bam 2>&1");
	}
	open COUNT, "<total.count" or die $!;
	my $header = <COUNT>;
	$header = <COUNT>;
	while(my $row = <COUNT>){
		chomp $row;
		my @cols = split /\t/, $row;
		$count{$cols[0]}{"exon"} = $cols[6];
		$count{$cols[0]}{"length"} = $cols[5];
		$countSum += $cols[6];
	}
	close COUNT;

	print $main::tee "\nRead Count:".$countSum."\n";

	open OUT, ">".$tag.".txt" or die $!;
	print OUT "Gene\tCount\tLength\n";
	foreach my $name (sort keys %count){
		print OUT "$name\t$count{$name}{exon}\t$count{$name}{length}\n";
	}
	close OUT;
}

1;

__END__
