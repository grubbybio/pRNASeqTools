#!/usr/bin/env perl
package Function;

use Modern::Perl;
use Function;
use File::Path qw(remove_tree make_path);

sub SRR { #download SRA files
	my ($self, $srr, $thread) = @_;
  if($srr =~ /([SE]RR\d+)$/ && !-e $srr){
		$srr = $1;

  	print $main::tee "Downloading...\n";

  	system ("fasterq-dump --threads ".$thread." --split-3 ".$srr);
  	if(-e $srr."_1.fastq"){
			unlink $srr.".fastq" if (-e $srr.".fastq");
  		return ($srr."_1.fastq", $srr."_2.fastq");
  	}else{
  		return ($srr.".fastq");
  	}
  }eles{
    return $srr;
  }
}

sub sum {
	my $self = shift @_;
	my $total = 0;
	foreach (@_) {
		$total += $_;
	}
	return $total;
}

sub average {
	my $self = shift @_;
	my $total = Function->sum(@_);
	my $average = $total / @_;
	return $average;
}

sub stdev {
	my $self = shift @_;
	my $average = Function->average(@_);
	my $sqtotal = 0;
	foreach(@_) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@_-1)) ** 0.5;
	return $std;
}

sub revcomp {
  my ($self, $inseq) = @_;
  my $revcomp = reverse($inseq);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub unzip { #seq file format check and unzip
  my ($self, $file, $tag) = @_;
  if(-e $file){
    if($file =~ /bz2$/){

      print $main::tee "Decompressing...\n";

      system ("bzip2 -dc ".$file." > ".$tag.".fastq");
    }elsif($file =~ /gz$/){

      print $main::tee "Decompressing...\n";

      system ("gzip -dc ".$file." > ".$tag.".fastq");
    }elsif($file =~ /fastq$/ or $file =~ /fq$/){
      if ($file ne $tag.".fastq"){

        print $main::tee "Renaming...\n";

        system ("cp ".$file." ".$tag.".fastq");
        unlink $file if($file =~ /^SRR/);
      }else{

        print $main::tee "Backing up...\n";

        system ("cp ".$file." ".$tag.".fastq.bak");
      }
    }else{
			die "Please provide the seq file in correct formats!"
		}

    print $main::tee "Completed!\n";

  }else{
    die "Please provide the correct seq file!";
  }
}

sub rmvc { #remove 3' polyC and reverse compliment
	my ($self, @tag) = @_;
	my $trig = 1;
	my ($clen, $seqname, $outseq, $outq, %rmvc, %rmvc2, $fiseq, $index);

	print $main::tee "Start to remove the 3' PolyC...\n";

	open TRIM, $tag[0]."_trimmed.fastq" or die $!;
	while (my $row = <TRIM>){
		chomp $row;
		if($. % 4 == 0){
			$rmvc{$seqname}{"outq"} = reverse substr($row, 0, -$clen) if(exists $rmvc{$seqname}{"outseq"});
		}elsif($. % 4 == 2){
			if($row =~ /^([A-Z]+[ATG])(C+)$/){
				$rmvc{$seqname}{"outseq"} = Function->revcomp($1);
				$clen = length ($2);
			}
		}elsif($. % 4 == 1){
			($seqname, $index) = split / /, $row;
		}
	}
	close TRIM;
	if($tag[1]){
		open TRIM, $tag[1]."_trimmed.fastq" or die $!;
		while (my $row = <TRIM>){
			chomp $row;
			if($. % 4 == 0){
				$rmvc2{$seqname}{"outq"} = reverse substr($row,$clen) if(exists $rmvc2{$seqname}{"outseq"});
			}elsif($. % 4 == 2){
				if($row =~ /^(G+)([ATC][A-Z]+)$/){
					$clen = length ($1);
					$rmvc2{$seqname}{"outseq"} = Function->revcomp($2);
				}
			}elsif($. % 4 == 1){
				($seqname, $index) = split / /, $row;
			}
		}
		close TRIM;
	}
	open COUT, ">".$tag[0].".fastq" or die $!;
	open COUT2, ">".$tag[1].".fastq" or die $! if($tag[1]);
	foreach $seqname (keys %rmvc){
		if($tag[1]){
			if(exists $rmvc{$seqname}{"outq"} and exists $rmvc2{$seqname}{"outq"}){
				print COUT "$seqname\n$rmvc{$seqname}{outseq}\n+\n$rmvc{$seqname}{outq}\n";
				print COUT2 "$seqname\n$rmvc2{$seqname}{outseq}\n+\n$rmvc2{$seqname}{outq}\n";
			}
		}else{
			if(exists $rmvc{$seqname}{"outq"}){
				print COUT "$seqname\n$rmvc{$seqname}{outseq}\n+\n$rmvc{$seqname}{outq}\n";
			}
		}
	}
	close COUT;
	close COUT2 if($tag[1]);

	print $main::tee "The 3' PolyC is removed!\n";
}

sub countBAM {
	my ($self, $tag, $thread, $seqStrategy, $prefix, $genome) = @_;

	system ("samtools index ".$tag.".bam");
	system ("bamCoverage -b ".$tag.".bam -bs 5 -p ".$thread." --filterRNAstrand forward --normalizeUsing RPKM -o ".$tag.".forward.bw");
	system ("bamCoverage -b ".$tag.".bam -bs 5 -p ".$thread." --filterRNAstrand reverse --normalizeUsing RPKM -o ".$tag.".reverse.bw");

	print $main::tee "\nStart counting...\n";

	my %count = ();
	my $countSum = 0;
	if($seqStrategy eq "single"){
		system ("featureCounts -T ".$thread." -C -O -G ".$prefix."/reference/".$genome."_chr_all.fasta -s 0 -a ".$genome."_genes.gtf -o total.count ".$tag.".bam 2>&1");
	}elsif($seqStrategy eq "paired"){
		system ("featureCounts -T ".$thread." -p -B -C -O -G ".$prefix."/reference/".$genome."_chr_all.fasta -s 0 -a ".$genome."_genes.gtf -o total.count ".$tag.".bam 2>&1");
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
