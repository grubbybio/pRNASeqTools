#!/usr/bin/env perl
package mode::ribometh;

use Modern::Perl;
use MooseX::App::Command;
extends qw(mode);
use File::Path qw(remove_tree make_path);
use File::Cat;
use Function;
use Ref;
use validate_options;
use input;
use Cwd 'abs_path';

command_short_description q[Analysis for RiboMeth-seq];
command_long_description q[Analysis for RiboMeth-seq];
command_usage q[pRNASeqTools ribometh [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT1]=[file1]+[file2] ... --treatment [TREATMENT2]=[file1]+[file2] ... ];

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
option 'reference' => (
  is => 'rw',
  isa => 'Str',
  default => 'genome',
  documentation => q[reference transcriptome fasta file],
);
option 'readlength' => (
  is => 'rw',
  isa => 'Num',
  default => '50',
  documentation => q[Raw read length],
);
option 'coverage' => (
  is => 'rw',
  isa => 'Num',
  default => '1000',
  documentation => q[Minimal coverage for ScoreMean calculation],
);
option 'adaptor2' => (
  is => 'rw',
  isa => 'Str',
  default => '1',
  documentation => q[Adaptor for the second read in the pair],
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
  my $reference = $options{'reference'};
  my $readLength = $options{'readlength'};
  my $coverage = $options{'coverage'};
  my $adaptor2 = $options{'adaptor2'};

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

    print $main::tee "\nBuilding STAR reference index ...\n";
    my (%geneSeq, $geneID, %geneinfo);
    if($reference eq "genome"){
      %geneinfo = Ref->geneinfo($prefix, $genome);
      open REF, "<transcripts.fa" or die $!;
      while(my $ref = <REF>){
        chomp $ref;
        if($ref =~ />(\w+)\.(\d+)/){
          $geneID = $1;
          $geneSeq{$geneID}{"transcriptID"} = $2;
        }else{
          $geneSeq{$geneID}{"seq"} .= $ref;
        }
      }
      close REF;
      open REF, ">reference.fa" or die $!;
      open GFF, ">reference.gff" or die $!;
      foreach $geneID (sort keys %geneSeq){
        print REF ">$geneID\n$geneSeq{$geneID}{seq}\n";
        print GFF "$geneID\treference\tgene\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\n";
        print GFF "$geneID\treference\tmRNA\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\.$geneSeq{$geneID}{transcriptID};Parent=$geneID\n";
        print GFF "$geneID\treference\texon\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID:exon:1;Parent=$geneID\.$geneSeq{$geneID}{transcriptID}\n";
      }
      close GFF;
      close REF;
      unlink "transcripts.fa";
    }else{
      if($reference =~ /^~\/(.+)/){
        $reference = $ENV{"HOME"}."/".$1;
      }else{
        $reference = abs_path "../".$reference;
      }
      symlink $reference.".fa", "reference.fa";
      open REF, "<reference.fa" or die $!;
      while(my $ref = <REF>){
        chomp $ref;
        if($ref =~ />(.+)/){
          $geneID = $1;
          $geneSeq{$geneID}{"transcriptID"} = 1;
        }else{
          $geneSeq{$geneID}{"seq"} .= $ref;
        }
      }
      open GFF, ">reference.gff" or die $!;
      foreach $geneID (sort keys %geneSeq){
        $geneinfo{$geneID} = length $geneSeq{$geneID}{"seq"};
        print GFF "$geneID\treference\tgene\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\n";
        print GFF "$geneID\treference\tmRNA\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\.$geneSeq{$geneID}{transcriptID};Parent=$geneID\n";
        print GFF "$geneID\treference\texon\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID:exon:$geneSeq{$geneID}{transcriptID};Parent=$geneID\.1\n";
      }
      close REF;
      close GFF;
    }
    remove_tree "Genome" if(-e "Genome");
    make_path "Genome";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeFastaFiles reference.fa --sjdbGTFfile reference.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000 --genomeSAindexNbases 5");

    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];

  		print $main::tee "\nMapping $tag...\n";

  		if($file !~ /,/){
  			my @files = Function->SRR($file, $thread);
        if($#files == 0){
        	Function->unzip($files[0], $tag);
        	if(defined $adaptor){

            print $main::tee "Trimming...\n";

						system ("cutadapt -j ".$thread." -m 15 --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
          	rename $tag."_trimmed.fastq", $tag.".fastq";
        	}
          system ("STAR --genomeDir Genome --seedSearchStartLmax 15 --outSAMtype SAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
          unlink ($tag.".fastq");
        }else{
          Function->unzip($files[0], $tag."_R1");
          Function->unzip($files[1], $tag."_R2");
          if(defined $adaptor){

            print $main::tee "Trimming...\n";

            system ("cutadapt -j ".$thread." -m 15 --trim-n -a ".$adaptor." -A ".$adaptor2." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
            rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
            rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
          }
          system ("STAR --genomeDir Genome --seedSearchStartLmax 15 --outSAMtype SAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          unlink ($tag."_R1.fastq", $tag."_R2.fastq");
        }
  		}else{
       	my ($file1, $file2) = split /,/, $file;
				Function->unzip($file1, $tag."_R1");
				Function->unzip($file2, $tag."_R2");
  		  if(defined $adaptor){
          system ("cutadapt -j ".$thread." -m 15 --trim-n -a ".$adaptor." -A ".$adaptor2." -o ".$tag."_R1_trimmed.fastq -p ".$tag."_R2_trimmed.fastq ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
          rename $tag."_R1_trimmed.fastq", $tag."_R1.fastq";
          rename $tag."_R2_trimmed.fastq", $tag."_R2.fastq";
        }
        system ("STAR --genomeDir Genome --seedSearchStartLmax 15 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag."_R1.fastq ".$tag."_R2.fastq 2>&1");
        unlink ($tag."_R1.fastq", $tag."_R2.fastq");
  		}
      system ("samtools view -h Aligned.sortedByCoord.out.bam > ".$tag.".sam");
      cat 'Log.final.out', \*STDOUT;
      my %outputGene = ();
      open SAM, "<".$tag.".sam" or die $!;
      while(my $sam = <SAM>){
        chomp $sam;
        next if($sam =~ /^@/);
        my @row = split /\t/, $sam;
        $outputGene{$row[2]} ++;
      }
      close SAM;
      open REF, ">".$tag.".fa" or die $!;
      open GFF, ">".$tag.".gff" or die $!;
      foreach my $geneID (sort keys %outputGene){
        if($outputGene{$geneID} >= $coverage * length($geneSeq{$geneID}{"seq"}) / 2){
          print REF ">$geneID\n$geneSeq{$geneID}{seq}\n";
          print GFF "$geneID\treference\tgene\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\n";
          print GFF "$geneID\treference\tmRNA\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID\.$geneSeq{$geneID}{transcriptID};Parent=$geneID\n";
          print GFF "$geneID\treference\texon\t1\t$geneinfo{$geneID}\t.\t+\t.\tID=$geneID:exon:1;Parent=$geneID\.$geneSeq{$geneID}{transcriptID}\n";
        }
      }
      close GFF;
      close REF;
      my %end = ();
      open SAM, "<".$tag.".sam" or die $!;
      open FIL, ">".$tag.".filtered.sam" or die $!;
      while(my $sam = <SAM>){
        chomp $sam;
        if($sam =~ /^@/){
          print FIL "$sam\n";
          next;
        }
        my @row = split /\t/, $sam;
        if($outputGene{$row[2]} >= $coverage * length($geneSeq{$row[2]}{"seq"}) / 2){
          if($row[5] =~ /^(\d+)M$/){
            my $maplength = $1;
            my $bin = sprintf "%b", $row[1];
            if($bin =~ /1[01][01][01][01][01][01]$/ or $bin == 0){
              $end{$row[2]}{$row[3]} ++;
              $end{$row[2]}{$row[3]+$maplength} ++ if($maplength < $readLength);
            }elsif($bin =~ /1[01][01][01][01][01][01][01]$/){
              $end{$row[2]}{$row[3]+$maplength} ++;
              $end{$row[2]}{$row[3]-1} ++ if($maplength < $readLength);
            }
            print FIL "$sam\n";
          }
        }
      }
      close SAM;
      close FIL;
      open WIG, ">".$tag.".ends.wig" or die $!;
      foreach my $chr (sort keys %end){
        print WIG "variableStep chrom=$chr\n";
        foreach my $site (sort {$a <=> $b} keys %{$end{$chr}}){
          print WIG "$site\t$end{$chr}{$site}\n";
        }
      }
      close WIG;
      system ("samtools view -Sb $tag.filtered.sam > $tag.filtered.bam");
      system ("samtools index $tag.filtered.bam");
      unlink "Aligned.sortedByCoord.out.bam", $tag.".sam", $tag.".filtered.sam";

      print $main::tee "\nCalculating RiboMeth Scores...\n";

      system ("Rscript --vanilla ".$prefix."/scripts/RNAmodR.R ".$tag." ".$coverage);
    }
    unlink "reference.fa", "reference.gff", (glob "*.bt2"), (glob "Log.*"), "SJ.out.tab";
    remove_tree "Genome";
  }
}

1;

__END__
