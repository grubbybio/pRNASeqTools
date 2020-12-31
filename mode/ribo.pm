#!/usr/bin/env perl
package mode::ribo;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use validate_options;
use input;
use File::Path qw/remove_tree/;
use Function;
use Ref;

command_short_description q[Analysis for Ribo-seq];
command_long_description q[Analysis for Ribo-seq];
command_usage q[pRNASeqTools ribo [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT1]=[file1]+[file2] ... --treatment [TREATMENT2]=[file1]+[file2] ... ];

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
option 'mmap' => (
  is => 'rw',
  isa => 'Str',
  default => 'u',
  documentation => q[method for assigning multiple mapped reads. Allowed: u, n, f, r],
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
  my $mmap = $options{'mmap'};

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

    die "Please specify the 3\' adaptor!" if($adaptor eq "");
    my %geneinfo = Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    system ("bowtie-build -q ".$genome.".fa ".$genome);
    system ("gffread -T ".$prefix."/reference/".$genome."_genes.gff -o ".$genome.".gtf");

    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];

      $file = Function->SRR($file, $thread);
      Function->unzip($file, $tag);
      if(defined $adaptor){

        print $main::tee "\nTrimming $tag...\n";

        system ("cutadapt -j ".$thread." -m 18 -M 42 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
        rename $tag."_trimmed.fastq", $tag.".fastq";
      }

      print $main::tee "\nStart mapping...\n";

      system ("ShortStack --outdir ShortStack_".$tag." --align_only --bowtie_m 10 --ranmax 10 --mmap ".$mmap." --mismatches 0 --bowtie_cores ".$thread." --nohp --readfile ".$tag.".fastq --genomefile ".$genome.".fa 2>&1");

      print $main::tee "\nAlignment Completed!\n";

      system ("samtools view -h ShortStack_".$tag."/".$tag."_trimmed.bam | awk '{if(\$10!=\"*\" && \$3!=\"*\") print}' > ".$tag.".sam");
      system ("samtools view -Sb --thread ".$thread." ".$tag.".sam > ".$tag.".bam");
      system ("samtools index ".$tag.".bam");
      system ("bamCoverage -b ".$tag.".bam -bs 5 -p ".$thread." --ignoreDuplicates --filterRNAstrand forward --normalizeUsing RPKM -o ".$tag.".forward.bw");
      system ("bamCoverage -b ".$tag.".bam -bs 5 -p ".$thread." --ignoreDuplicates --filterRNAstrand reverse --normalizeUsing RPKM -o ".$tag.".reverse.bw");
      remove_tree("ShortStack_".$tag);
      unlink ($tag.".sam", $tag.".fastq");
    }

    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par) if(!$mappingonly);
    unlink ($genome.".fa", glob ($genome."*.ebwt"), $genome.".gtf", $genome.".fa.fai");

  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".bam", $pre.".bam" or die $!;
    }
    Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    system ("gffread -T ".$prefix."/reference/".$genome."_genes.gff -o ".$genome.".gtf");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par);
    unlink ($genome.".gtf", glob ("*.bam"), $genome.".fa");
  }
}

1;

__END__
