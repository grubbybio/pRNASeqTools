#!/usr/bin/env perl
package mode::degradome;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Path qw/remove_tree make_path/;
use File::Copy::Recursive qw/fmove/;
use File::Cat;
use Function;
use Ref;
use validate_options;
use input;
use mode::degradome;

command_short_description q[Peak finder for degradome sequencing (GMUCT or PARE)];
command_long_description q[Peak finder for degradome sequencing (GMUCT or PARE)];
command_usage q[pRNASeqTools degradome [OPTIONS] --control [CONTROL]=[file1]+[file2] ... ];

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
option 'targets' => (
  is => 'rw',
  isa => 'Str',
  default => 'all',
  documentation => q[Transcript list of interested targets for CRI. If not defined, all transcipts will be used],
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
  my $targets = $options{'targets'};

  my ($tags_ref, $files_ref, $par_ref) = input->run($control);
  my @tags = @$tags_ref;
  my @files = @$files_ref;
  if (defined $treatment){
    ($tags_ref, $files_ref, $par_ref) = input->run($treatment);
    push @tags, @$tags_ref;
    push @files, @$files_ref;
  }
  my $par = join " ", @tags;

  if(!$nomapping){

    Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    system ("gffread -T ".$prefix."/reference/".$genome."_genes.gff -o ".$genome.".gtf");
    remove_tree "Genome" if(-e "Genome");
    make_path "Genome";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeFastaFiles ".$genome.".fa --limitGenomeGenerateRAM 64000000000");
    remove_tree "Genome2" if(-e "Genome2");
    make_path "Genome2";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome2 --runMode genomeGenerate --genomeFastaFiles ".$prefix."/reference/".$genome."_chr_all.fasta --sjdbGTFfile ".$prefix."/reference/".$genome."_genes.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000");
    for(my $i=0;$i<=$#tags;$i++){
      my $file = $files[$i];
      my $tag = $tags[$i];

    	print $main::tee "\nWorking on $tag...\n";

      $file = Function->SRR($file, $thread);
      Function->unzip($file, $tag);
      if(defined $adaptor){

        print $main::tee "\nTrimming...\n";

        system ("cutadapt -j ".$thread." -m 18 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
        rename $tag."_trimmed.fastq", $tag.".fastq";
      }

      print $main::tee "\nStart mapping...\n";

      system ("STAR --genomeDir Genome --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
      rename "Aligned.sortedByCoord.out.bam", $tag.".bam";
      cat 'Log.final.out', \*STDOUT;
      system ("samtools index ".$tag.".bam");
      system ("bamCoverage -b ".$tag.".bam --skipNAs -bs 5 -p ".$thread." --minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o ".$tag.".bw");

      system ("STAR --genomeDir Genome2 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
      rename "Aligned.sortedByCoord.out.bam", $tag."genomic.bam";
      cat 'Log.final.out', \*STDOUT;
      system ("samtools index ".$tag."genomic.bam");
      system ("bamCoverage -b ".$tag."genomic.bam --skipNAs -bs 5 -p ".$thread." --minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o ".$tag.".genomic.bw");

      print $main::tee "\nAlignment Completed!\n";

      open FQ, "<".$tag.".fastq" or die $!;
      my %fqDemultiplex;
      while(my $fq = <FQ>){
        chomp $fq;
        if($. % 4 == 2){
          $fqDemultiplex{$fq} += 1;
        }
      }
      close FQ;
      open LIB, ">".$tag."_lib.txt" or die $!;
      foreach my $fqSeq (keys %fqDemultiplex){
        print LIB "$fqSeq\t$fqDemultiplex{$fqSeq}\n";
      }
      close LIB;
      unlink $tag.".fastq";
    }
    if(!$mappingonly){

      print $main::tee "Finding peaks...\n";

      mode::degradome->Peak($thread, $prefix, $genome, \@tags);

      print $main::tee "Calculating CRIs...\n";

      system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par);
      mode::degradome->CRI(@tags, $targets);
      system ("Rscript --vanilla ".$prefix."/scripts/CRI.R ".$par);
    }
    unlink glob ("Log.*"), "SJ.out.tab", $genome.".fa", $genome.".gtf", $genome.".fa.fai";
    remove_tree ("Genome", "Genome2");
  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".bam", $pre.".bam" or die $!;
      symlink "../".$pre."_lib.txt", $pre."_lib.txt";
    }

    print $main::tee "Finding peaks...\n";

    Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    mode::degradome->Peak($thread, $prefix, $genome, \@tags);

    print $main::tee "Calculating CRIs...\n";

    system ("gffread -T ".$prefix."/reference/".$genome."_genes.gff -o ".$genome.".gtf");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par);
    mode::degradome->CRI(@tags, $targets);
    system ("Rscript --vanilla ".$prefix."/scripts/CRI.R ".$par);
    unlink (glob ("*_lib.txt"), $genome.".gtf", glob ("*.bam"), $genome.".fa");
  }
}

sub CRI {
	my ($self, @tags) = @_;
	my $targets = pop @tags;

	for(my $i=0;$i<=$#tags;$i++){
		my $tag = $tags[$i];
		open RIBO, "<".$tag.".csv" or die $!;
		my $dummy = <RIBO>;
		my %frame = ();
		while(my $ribo = <RIBO>){
			chomp $ribo;
			my @read = split /,/, $ribo;
			if($read[9] eq "cds"){
				my $tmp2 = $read[1] - $read[5];
				if($tmp2 > 0){
					$frame{$read[0]}{$tmp2 % 3} ++;
					$frame{$read[0]}{"total"} ++;
				}
			}
		}
		close RIBO;
		open CRI, ">".$tag."_CRI.txt" or die $!;
		if($targets eq "all"){
			foreach my $transcript (sort keys %frame){
				if($frame{$transcript}{"total"} >= 100){
          $frame{$transcript}{0} = 0 if(!exists $frame{$transcript}{0});
          $frame{$transcript}{1} = 0 if(!exists $frame{$transcript}{1});
          $frame{$transcript}{2} = 0 if(!exists $frame{$transcript}{2});
					my $log = log(($frame{$transcript}{1}+$frame{$transcript}{2}+1)/(2*$frame{$transcript}{0}+1))/log(2);
					print CRI "$transcript\t$frame{$transcript}{0}\t$frame{$transcript}{1}\t$frame{$transcript}{2}\t$log\n";
				}
			}
		}else{
			open TAR, "<".$targets or die $!;
			while(my $target = <TAR>){
				if(exists $frame{$target}{"total"} and $frame{$target}{"total"} >= 100){
					my $log = log(($frame{$target}{1}+$frame{$target}{2}+1)/(2*$frame{$target}{0}+1))/log(2);
					print CRI "$target\t$frame{$target}{0}\t$frame{$target}{1}\t$frame{$target}{2}\t$log\n";
				}
			}
			close TAR;
		}
		close CRI;
	}
}

sub Peak {
	my ($self, $thread, $prefix, $genome, $tags_ref) = @_;
  my @tags = map {$_."_lib.txt"} @$tags_ref;
  my $tags_lib = join " ", @tags;
  make_path "sparta";
  chdir "sparta";
  symlink "../".$genome.".fa", $genome.".fa";
  my %mirna;
  my %fas = Ref->fas($prefix, $genome);
  open MI, "<".$prefix."/reference/".$genome."_miRNA_miRNA_star.gff" or die $!;
  while (my $mi = <MI>){
    chomp $mi;
    my @row = split /\t/, $mi;
    my $miseq = substr($fas{$row[0]}, $row[3]-1, $row[4]-$row[3]+1);
    $miseq = Function->revcomp($miseq) if($row[6] eq "-");
    $mirna{$miseq}{"name"} .= $row[8].";";
  }
  close MI;
  open MIR, ">".$genome."_miRNA.fa" or die $!;
  foreach my $miseq (sort keys %mirna){
    chop $mirna{$miseq}{"name"};
    print MIR ">$mirna{$miseq}{name}\n$miseq\n";
  }
  close MIR;
  symlink $prefix."/sPARTA.py", "sPARTA.py";
  foreach my $file (@tags){
    symlink "../".$file, $file or die $!;
  }
  system ("python3 sPARTA.py -accel ".$thread." -featureFile ".$genome.".fa -genomeFeature 0 -miRNAFile ".$genome."_miRNA.fa -libs ".$tags_lib." -minTagLen 18 -tarPred -tarScore --tag2FASTA --map2DD --validate");
  opendir my $output, "output" or die $!;
  foreach my $file (readdir $output){
    fmove "output/".$file, "..";
  }
  chdir "..";
  remove_tree "sparta";
}

1;

__END__
