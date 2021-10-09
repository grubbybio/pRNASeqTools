#!/usr/bin/env perl
package mode::ribo;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use validate_options;
use input;
use File::Path qw/remove_tree make_path/;
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
option 'mismatches' => (
  is => 'rw',
  isa => 'Num',
  default => '1',
  documentation => q[Number of mismatches allowed. Allowed: 0, 1, 2],
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
  my $mismatches = $options{'mismatches'};

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
    remove_tree "Genome" if(-e "Genome");
    make_path "Genome";
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeSAindexNbases 12 --genomeFastaFiles ".$prefix."/reference/".$genome."_chr_all.fasta --sjdbGTFfile ".$prefix."/reference/".$genome."_genes.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000");
    system ("gffread -T -C -o ".$genome.".gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    system ("grep \"^chr[0-9]\" ".$genome.".gtf > tmp; mv tmp ".$genome.".gtf");
    for(my $i=0;$i<=$#tags;$i++){
      my $tag = $tags[$i];
      my $file = $files[$i];
      $file = Function->SRR($file, $thread);
      Function->unzip($file, $tag);
      if(defined $adaptor){
        print $main::tee "\nTrimming $tag...\n";
        system ("cutadapt -j ".$thread." -m 18 --discard-untrimmed --trim-n -a ".$adaptor." -o ".$tag."_trimmed.fastq ".$tag.".fastq 2>&1");
        rename $tag."_trimmed.fastq", $tag.".fastq";
      }
      print $main::tee "\nStart mapping...\n";
      system ("STAR --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --quantMode TranscriptomeSAM --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
      system ("samtools view -h -F 0x100 Aligned.toTranscriptome.out.bam |samtools sort -o ".$tag.".bam -");
      system ("samtools index ".$tag.".bam");
      system ("bamCoverage -b ".$tag.".bam -bs 5 -p ".$thread." --filterRNAstrand forward --normalizeUsing RPKM -o ".$tag.".bw");
      unlink $tag.".fastq";
    }
    unlink ("Aligned.sortedByCoord.out.bam", "Aligned.toTranscriptome.out.bam", glob ("Log.*"), "SJ.out.tab");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par) if(!$mappingonly);
    unlink ($genome.".gtf");
    remove_tree "Genome";
  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".bam", $pre.".bam" or die $!;
    }
    Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    system ("gffread -T -C -o ".$genome.".gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par);
    unlink ($genome.".gtf", glob ("*.bam"), $genome.".fa");
  }
}

1;

__END__
