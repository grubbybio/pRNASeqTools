#!/usr/bin/env perl
package mode::clip;

use Modern::Perl;
use MooseX::App::Command;
extends qw(mode);
use Function;
use File::Path qw(remove_tree make_path);
use File::Cat;
use validate_options;
use input;
use mode::clip;

command_short_description q[Analysis for CLIP-seq];
command_long_description q[Analysis for CLIP-seq];
command_usage q[pRNASeqTools clip [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT]=[file1]+[file2] ... ];

option 'no-mapping' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Just perform the statistic analysis.],
);
option 'mapping-only' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Do not perform the statistic analysis.],
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
  default => '0.05',
  documentation => q[Threshold for DEG in P value.],
);

sub run {
  my ($self) = @_;

  my %options = validate_options->run($self);
  my $thread = $options{'thread'};
  my $genome = $options{'genome'};
  my $adaptor = $options{'adaptor'};
  my $prefix = $options{'prefix'};
  my $control = $options{'control'};
  my $treatment = $options{'treatment'};
  my $nomapping = $options{'no-mapping'};
  my $mappingonly = $options{'mapping-only'};
  my $foldchange = $options{'foldchange'};
  my $pvalue = $options{'pvalue'};

  my ($tags_ref, $files_ref, $par_ref) = input->run($control);
  my @tags = @$tags_ref;
  my @files = @$files_ref;
  my @par = @$par_ref;
  if (defined $treatment){
    my ($tags_ref, $files_ref, $par_ref) = input->run($treatment);
    push @tags, @$tags_ref;
    push @files, @$files_ref;
    push @par, @$par_ref;
  }

  if(!$nomapping){
		die "Please specify the 3\' adaptor!" if($adaptor eq "");

    print $main::tee "\nBuilding STAR genome index ...\n";

    remove_tree "Genome" if(-e "Genome");
    make_path "Genome";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeFastaFiles ".$prefix."/reference/".$genome."_chr_all.fasta --sjdbGTFfile ".$prefix."/reference/".$genome."_genes.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID");

    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];

      print $main::tee "\nMapping $tag...\n";

			if($file !~ /,/){
				my @files = Function->SRR($file, $thread);
	      if($#files == 0){
					Function->unzip($files[0], $tag);
					system ("cutadapt -j ".$thread." -m 18 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
		  		system ("STAR --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outReadsUnmapped Fastx --outSAMmultNmax 1 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_trimmed.fastq 2>&1");
		  		rename "Unmapped.out.mate1", $tag.".unmapped.fastq";
          unlink ($tag.".fastq", $tag."_trimmed.fastq");
	      }else{
	        Function->unzip($files[0], $tag."_R1");
	        Function->unzip($files[1], $tag."_R2");
	        system ("cutadapt -j ".$thread." -m 18 --trim-n -a ".$adaptor." -A AAAAAAGAAAAAA -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
	        system ("STAR --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outReadsUnmapped Fastx --outSAMmultNmax 1 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1_trimmed.fastq ".$tag."_R2_trimmed.fastq 2>&1");
	        rename "Unmapped.out.mate1", $tag.".unmapped_R1.fastq";
          rename "Unmapped.out.mate2", $tag.".unmapped_R2.fastq";
          unlink ($tag."_R1.fastq", $tag."_R2.fastq", $tag."_R1_trimmed.fastq", $tag."_R2_trimmed.fastq");
	      }
			}else{
	     	my ($file1, $file2) = split /,/, $file;
				Function->unzip($file1, $tag."_R1");
				Function->unzip($file2, $tag."_R2");
	      system ("cutadapt -j ".$thread." -m 18 --trim-n -a ".$adaptor." -A AAAAAAGAAAAAA -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
	      system ("STAR --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outReadsUnmapped Fastx --outSAMmultNmax 1 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1_trimmed.fastq ".$tag."_R2_trimmed.fastq 2>&1");
				rename "Unmapped.out.mate1", $tag.".unmapped_R1.fastq";
        rename "Unmapped.out.mate2", $tag.".unmapped_R2.fastq";
        unlink ($tag."_R1.fastq", $tag."_R2.fastq", $tag."_R1_trimmed.fastq", $tag."_R2_trimmed.fastq");
			}
			rename "Aligned.sortedByCoord.out.bam", $tag.".bam";
      cat 'Log.final.out', \*STDOUT;
      system ("samtools index ".$tag.".bam");
      system ("bamCoverage -b ".$tag.".bam --skipNAs -bs 5 -p ".$thread." --ignoreDuplicates --normalizeUsing CPM -o ".$tag.".bw");

      print $main::tee "\nAlignment completed!\nFinding peaks...\n";

      system ("clipper -b ".$tag.".bam -s ".$genome." --FDR=0.01 --minreads=2 --processors=".$thread." --threshold-method=binomial --min_width=20 -o ".$tag.".fitted_clusters.bed -v 2>&1");

    }
    unlink ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab");
    remove_tree "Genome";
    if(!$mappingonly && $#par > 1){
      mode::clip->ana(\@par, $prefix, $pvalue, $foldchange);
    }
  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".nf", $pre.".nf" or die $!;
      symlink "../".$pre.".bam", $pre.".bam";
      symlink "../".$pre.".fitted_clusters.bed", $pre.".fitted_clusters.bed";
    }

    mode::clip->ana(\@par, $prefix, $pvalue, $foldchange);
    unlink glob ("*.bed"), glob ("*.nf"), glob ("*.bam*");
  }
}

sub ana {
  my ($self, $pars_ref, $prefix, $pvalue, $foldchange) = @_;

  print $main::tee "\nOverlaping peaks...\n";

  my @pars = @$pars_ref;
  my $par = join " ", @pars;
  my (@tags, $command);
  while (my @group = splice @pars, 0, 2){
    for(my $s=1;$s<=$group[1];$s++){
      push @tags, $group[0]."_".$s;
      $command .= $group[0]."_".$s.".fitted_clusters.bed ";
    }
  }
  system ("cat ".$command."\| sort -k1,1 -k2,2n \| bedtools merge -c 4,6 -o collapse,distinct -s -i - > tmp.bed");
  open PEAK, "<tmp.bed" or die $!;
  open FI, ">ref.bed" or die $!;
  while(my $dd = <PEAK>){
    chomp $dd;
    my $sum = 0;
    my @sample = ();
    my %PGene = ();
    my @row = split /\t/, $dd;
    my @peaks = ();
    if($row[3] =~ /,/){
      @peaks = split /,/, $row[3];
    }else{
      @peaks = ($row[3]);
    }
    foreach my $peak (@peaks){
      @sample = split /_/, $peak;
      $PGene{$sample[0]} = "";
      $sum += $sample[2];
    }
    my $OutputGene = join "_", sort keys %PGene;
    print FI "$row[0]\t$row[1]\t$row[2]\t$OutputGene\t$row[3]\_$sum\t$row[4]\n" if $sum >= 10;
  }
  close PEAK;
  close FI;

  print $main::tee "Counting reads in each sample...\n";

  foreach my $tag (@tags){
    system ("bedtools intersect -c -s -a ref.bed -b ".$tag.".bam |awk \'{print \$1\"_\"\$2\"_\"\$3\"_\"\$6\"_\"\$4\"\\t\"\$7}\' > ".$tag.".txt");
  }

  print $main::tee "Finding differential peaks...\n";

  system ("Rscript --vanilla ".$prefix."/scripts/CLIP.R ".$pvalue." ".$foldchange." ".$par);
  unlink "tmp.bed", "ref.bed";
}

1;

__END__
