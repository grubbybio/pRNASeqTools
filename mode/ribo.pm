#!/usr/bin/env perl
package mode::ribo;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Cat;
use validate_options;
use input;
use File::Path qw/remove_tree make_path/;
use Function;
use Ref;
use File::Copy qw/mv/;
use Cwd qw/abs_path/;

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
  my $mask = $options{'mask'};

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
  	if(defined $mask){
  		if($mask =~ /^~\/(.+)/){
  			$mask = $ENV{"HOME"}."/".$1;
  		}elsif($mask !~ /^\//){
  			$mask = abs_path "../".$mask;
  		}
  		symlink $mask, "mask.fa";
  		system ("bowtie-build -q mask.fa mask");
  	}
    remove_tree "Genome" if(-e "Genome");
    make_path ("Genome", "chr");
    system ("gffread -T -C -o ".$genome.".gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    Ref->PrimaryTranscript($prefix, $genome);
    system ("STAR --runThreadN ".$thread." --genomeDir Genome --runMode genomeGenerate --genomeSAindexNbases 12 --genomeFastaFiles ".$prefix."/reference/".$genome."_chr_all.fasta --sjdbGTFfile ".$genome.".PrimaryTranscript.gtf --limitGenomeGenerateRAM 64000000000");
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
      if(defined $mask){
        system ("bowtie -v 2 -a --un tmp.fastq -p ".$thread." -t -x mask ".$tag.".fastq ".$tag.".mask.out 2>&1");
        rename "tmp.fastq", $tag.".fastq";
        unlink $tag.".mask.out";
      }
      print $main::tee "\nStart mapping...\n";
      system ("STAR --genomeDir Genome --alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outSAMmultNmax 1 --outFilterMultimapNmax 50 --quantMode TranscriptomeSAM --outFilterMismatchNoverLmax 0.1 --runThreadN ".$thread." --readFilesIn ".$tag.".fastq 2>&1");
      cat 'Log.final.out', \*STDOUT;
      system ("samtools view -h -F 0x100 Aligned.toTranscriptome.out.bam |samtools sort -o ".$tag.".bam -");
      system ("samtools index ".$tag.".bam");
      system ("samtools view -q 10 -b Aligned.sortedByCoord.out.bam > ".$tag.".chr.bam");
      system ("samtools index ".$tag.".chr.bam");
      system ("bamCoverage -b ".$tag.".chr.bam --skipNAs -bs 5 --minMappingQuality 10 -p ".$thread." --normalizeUsing CPM -o ".$tag.".bw");
      unlink ($tag.".fastq", "tmp.sam");
      mv $tag.".chr.bam", "chr";
      mv $tag.".chr.bam.bai", "chr";
      mv $tag.".bw", "chr";
    }
    unlink ("Aligned.toTranscriptome.out.bam", "Aligned.sortedByCoord.out.bam", glob ("Log.*"), "SJ.out.tab", $genome.".gtf");
    unlink ("mask.fa", glob ("mask*ebwt")) if(-e "mask.fa");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par) if(!$mappingonly);
    remove_tree "Genome";
  }else{
    foreach my $pre (@tags){
      symlink "../".$pre.".bam", $pre.".bam" or die $!;
    }
    Ref->geneinfo($prefix, $genome);
    rename "transcripts.fa", $genome.".fa";
    Ref->PrimaryTranscript($prefix, $genome);
    system ("gffread -T -C -o ".$genome.".gtf -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
    system ("Rscript --vanilla ".$prefix."/scripts/ribo.R ".$genome." ".$par);
    unlink ($genome.".gtf", glob ("*.bam"), $genome.".fa");
  }
}

1;

__END__
