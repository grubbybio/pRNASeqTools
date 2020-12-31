#!/usr/bin/env perl
package mode::phasi;

use Modern::Perl;
use MooseX::App::Command;
extends qw/mode/;
use File::Path qw/remove_tree/;
use Function;
use Ref;
use validate_options;
use input;

command_short_description q[Analysis for phasiRNA using small RNA-seq];
command_long_description q[Analysis for phasiRNA using small RNA-seq];
command_usage q[pRNASeqTools phasi [OPTIONS] --control [CONTROL]=[file1]+[file2] ... --treatment [TREATMENT]=[file1]+[file2] ... ];

option 'no-mapping' => (
  is => 'rw',
  isa => 'Bool',
  default => 0,
  documentation => q[Just perform the statistic analysis],
);
option 'mmap' => (
  is => 'rw',
  isa => 'Str',
  default => 'u',
  documentation => q[method for assigning multiple mapped reads. Allowed: u, n, f, r],
);
option 'norm' => (
  is => 'rw',
  isa => 'Str',
  default => 'rRNA,total',
  documentation => q[method for normalization, seperated by comma. Allowed: rRNA, total.],
);
option 'period' => (
	is => 'rw',
	isa => 'Num',
	default => 21,
	documentation => q[period size. Allowed: 19 - 26],
);
option 'binsize' => (
	is => 'rw',
  isa => 'Num',
  default => 100,
  documentation => q[windows size to slice the genome.],
);
option 'phasingscore' => (
	is => 'rw',
	isa => 'Num',
	default => 50,
	documentation => q[phasing score cutoff for qualified bins.],
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
  my $binsize = $options{'binsize'};
  my $norm = $options{'norm'};
  my @norms = split /,/, $norm;
  my $period = $options{'period'};
  my $phasingscore = $options{'phasingscore'};

  my ($tags_ref, $files_ref, $par_ref) = input->run($control);
  my @tags = @$tags_ref;
  my @files = @$files_ref;
  my @pars = @$par_ref;
  if (defined $treatment){
    ($tags_ref, $files_ref, $par_ref) = input->run($treatment);
    push @tags, @$tags_ref;
    push @files, @$files_ref;
    push @pars, @$par_ref;
  }

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
      system ("bowtie -v 2 -a -p ".$thread." -t ".$prefix."/reference/lsu_rrna ".$tag.".fastq ".$tag.".rRNA.out 2>&1");
      system ("awk -F \"\t\" \'BEGIN{x=0}{if(/$genome/){x++}}END{print \"rRNA\\t\"x}\' ".$tag.".rRNA.out > ".$tag.".nf");

      system ("ShortStack --outdir ShortStack_".$tag." --align_only --bowtie_m 1000 --ranmax 50 --mmap ".$mmap." --mismatches 0 --bowtie_cores ".$thread." --nohp --readfile ".$tag.".fastq --genomefile ".$prefix."/reference/".$genome."_chr_all.fasta 2>&1");
      unlink (glob ($tag."*.fastq"), $tag.".rRNA.out");

      print $main::tee "\nAlignment Completed!\n";

      system ("samtools view -h ShortStack_".$tag."/".$tag.".bam > ".$tag);
      system ("awk '{if(\$0~/^@/) print > (FILENAME\".unmapped.sam\"); if(\$10!=\"*\" && \$3!=\"*\") print > (FILENAME\".sam\"); if(\$10!=\"*\" && \$3==\"*\") print > (FILENAME\".unmapped.sam\")}' ".$tag);
      system ("samtools view -Sb ".$tag.".unmapped.sam > ".$tag.".unmapped.bam");
      system ("samtools view -Sb ".$tag.".sam > ".$tag.".bam");
      system ("awk \'{n++}END{print \"total\\t\"n}\' ".$tag.".sam >> ".$tag.".nf");
      unlink ($tag."", $tag.".unmapped.sam", $tag.".sam");
      remove_tree "ShortStack_".$tag;
    }
    mode::phasi->phasi(\@pars, $period, \@norms, $prefix, $genome, $binsize, $phasingscore);
  }else{
    foreach my $pre (@tags){
    	symlink "../".$pre.".nf", $pre.".nf" or die $!;
    	symlink "../".$pre.".bam", $pre.".bam";
    }
    mode::phasi->phasi(\@pars, $period, \@norms, $prefix, $genome, $binsize, $phasingscore);
    unlink glob ("*.bam"), glob ("*.nf");
  }
}

sub phasi {
  my ($self, $par_ref, $period, $norm_ref, $prefix, $genome, $binsize, $phasingscore) = @_;
  my @par = @$par_ref;
  my @norms = @$norm_ref;
  while (my @group = splice @par, 0, 2){
  	my $sample = $group[0];
  	my $num = $group[1];

		print $main::tee "\nphasiRNA analysis...\n";

		my $count = 0;
  	my (%data1, %data2, $rc, %gene_len, %phase_res, %hash2, %normc, %data0);
  	my $command = "samtools merge ".$sample.".bam ";
  	for(my $j=1;$j<=$num;$j++){
  		open my $rrna, "<$sample\_$j.nf" or die $!;
	    while (my $r = <$rrna>){
	      chomp $r;
	      my ($rr, $base) = split /\t/, $r;
	      $normc{$rr} += $base;
	    }
    	close $rrna;
    	$command .= $sample."_".$j.".bam ";
	  }
	  system $command;
	  system ("samtools view -h ".$sample.".bam > ".$sample.".sam");

	  print $main::tee "Reading SAM ...\n";

	  open my $sam, "$sample.sam" or die $!;
	  open PHASAM, ">$sample.phasi.sam" or die $!;
	  while (my $ee = <$sam>){
	    chomp $ee;
	    if ($ee =~ /^\s*$/){
	      print PHASAM "$ee\n";
	      next;
	    }
	    if ($ee =~ /^\@/){
	      print PHASAM "$ee\n";
	      next;
	    }
	    my ($readID, $flag, $chr, $pos, $MAPQ, $len, $line) = split /\t/, $ee;
	    $data0{$chr}{$pos} ++;
	    if($len eq $period."M"){
		    print PHASAM "$ee\n";
		    if ($flag == 0){
		      $data1{$chr}{$pos}{0} ++;
		      $count ++;
		    }  #forward strand
		    elsif ($flag == 16) {
		      $data2{$chr}{$pos}{16} ++;
		      $count ++;
		    } #reverse strand
		  }
	    if($count > 0 && $count % 100000 == 0){
	      print STDERR "Read counts: $count\r";
	    }
	  }
	  close PHASAM;
	  system ("samtools view -Sb "."$sample".".phasi.sam > "."$sample".".phasi.bam");
	  unlink ($sample.".phasi.sam", $sample.".sam", $sample.".bam");
	  close $sam;

	  print $main::tee "Read counts: $count\nCombining strands ...\n";

	  for my $chr (keys %data2) {
	    for my $pos (keys %{$data2{$chr}}) {
	      # 2-nt overhang
	      if (exists $data1{$chr}{$pos + 2}{0}) {
	        $data1{$chr}{$pos+2}{0} += $data2{$chr}{$pos}{16};
	      } else {
	        $data1{$chr}{$pos+2}{0} = $data2{$chr}{$pos}{16};
	      }
	    }
	  }

  	foreach my $mnorm (@norms){
  		my %data3 = ();
			my $rc = $normc{$mnorm};
		  for my $chr (keys %data1){
		    for my $pos (keys %{$data1{$chr}}) {
		        $data3{$chr}{$pos}{0} = $data1{$chr}{$pos}{0} * 1000000 / $rc;
		    }
		  }

		  print $main::tee "Calculating phasing score ...\n";

		  for my $chr (sort keys %data3) {

		    print $main::tee "$chr ...\n";

		    for my $pos (sort {$a <=> $b} keys %{$data3{$chr}}) {
		      my ($n, $phased, $unphased, $phasing_score, $xx) = (0, 0, 0, 0, 0);
		      for my $cycle (0..9) { # 10 cycles B.Meyers lab, RNA, 2009
		        my $new_pos = $pos + $period * $cycle; # phased position
		        $xx += $data0{$chr}{$new_pos} if(exists $data0{$chr}{$new_pos});
		        if (exists $data3{$chr}{$new_pos}{0}) {
		          $n ++;  # num of pos plus 1
		          $phased += $data3{$chr}{$new_pos}{0}; # phased cumulation
		        }
		        for my $unphased_nt_per_cycle (1..$period-1) {
		          $new_pos ++; # unphased positions
		          $xx += $data0{$chr}{$new_pos} if(exists $data0{$chr}{$new_pos});
		          if (exists $data3{$chr}{$new_pos}{0}) {
		            $n ++; # num of pos plus 1
		            $unphased += $data3{$chr}{$new_pos}{0};# unphased cumulation
		          }
		        }
		      }
		      my $ratio = ($phased + $unphased) / (($xx + 1) * 1000000 / $rc);
		      if ($n < 3 || $phased < 10 || $ratio < 0.25) {
		        $phasing_score = 0;
		      } else {
		        $phasing_score = ($n - 2) * log (1 + 10 * $phased / ($unphased + 1)); # B. Meyers lab 2009
		      }
		      $phase_res{$chr}{$pos} = $phasing_score;
		    }
		  }

		  print $main::tee "Outputing result ...\n";

		  open my $phase_res_out, ">$sample.$mnorm.phasiRNA.txt";
		  open my $bg, ">$sample.$mnorm.phasiRNA.bedgraph";
		  open my $plotbg, ">$sample.$mnorm.plot.phasiRNA.bedgraph";
		  print $phase_res_out "#BIN,Phase_score\n";
		  for my $chr (sort keys %phase_res){
		    for my $pos (sort {$a <=> $b} keys %{$phase_res{$chr}}){
		      if ($phase_res{$chr}{$pos} > 0){
		        print $plotbg "$chr\t$pos\t$pos\t$phase_res{$chr}{$pos}\n";
  		      if($phase_res{$chr}{$pos} >= $phasingscore and $data3{$chr}{$pos}{0} >= 1){
  		        print $bg "$chr\t$pos\t$pos\t$phase_res{$chr}{$pos}\n";
  		        my $bin = int($pos/$binsize);
  		        if(!exists $hash2{$chr}{$bin}){
  		          $hash2{$chr}{$bin} = $phase_res{$chr}{$pos};
  		        }elsif($phase_res{$chr}{$pos} > $hash2{$chr}{$bin}){
  		          $hash2{$chr}{$bin} = $phase_res{$chr}{$pos};
  		        }
  		      }
          }
		    }
		  }
		  for my $chr (sort keys %hash2){
		    for my $bin (sort {$a <=> $b} keys %{$hash2{$chr}}){
		      print $phase_res_out "$chr\_$bin,$hash2{$chr}{$bin}\n";
		    }
		  }
      %phase_res = ();
      %hash2 = ();
		  close $phase_res_out;
		  close $bg;
		  close $plotbg;

		  my %ann = Ref->ann($prefix, $genome, $binsize);
		  open TMP2, "$sample.$mnorm.phasiRNA.txt" or die $!;
		  open TMP3, ">tmp3" or die $!;
		  while(my $bb = <TMP2>){
		    chomp $bb;
		    my @row = split /,/, $bb;
		    $row[0] =~ s/"//g;
		    if(exists $ann{$row[0]}){
		      print TMP3 "$bb,$ann{$row[0]}\n";
		    }else{
		      print TMP3 "$bb,NA\n";
		    }
		  }
		  close TMP2;
		  close TMP3;
		  rename "tmp3", $sample.".".$mnorm.".phasiRNA.txt";

		  print $main::tee "\nphasiRNA analysis complete!\n\n";
		}
	}
}

1;

__END__
