#!/usr/bin/env perl
package Ref;

use Modern::Perl;
use Ref;
use Function;

sub gff {
  my ($self, $prefix, $genome) = @_;
  my $id = "";
  my ($idt, $ind, $chr, %index);
  open GFF, "$prefix/reference/".$genome."_genes.gff" or die $!;
  while(my $aa = <GFF>){
    chomp $aa;
    my @row = split /\t/, $aa;
    next if($row[2] =~ /UTR/ || $row[2] =~ /c_transcript/ || $row[2] =~ /region/);
    next if($row[2] eq "protein" || $row[2] eq "CDS" ||$row[2] =~ /[^i]RNA/);
    if($row[2] =~ /gene/){
      if($row[8] =~ /^ID=(\w+);/o){
        $id = $1;
        $ind = int($row[3]/100000);
        $chr = $row[0];
        $index{$chr}{$ind}{$id}{start} = $row[3];
        $index{$chr}{$ind}{$id}{end} = $row[4];
        $index{$chr}{$ind+1}{$id}{start} = $row[3];
        $index{$chr}{$ind+1}{$id}{end} = $row[4];
        $index{$chr}{$ind-1}{$id}{start} = $row[3];
        $index{$chr}{$ind-1}{$id}{end} = $row[4];
      }
    }elsif($row[8] =~ /Parent=$id/){
      $index{$chr}{$ind}{$id}{exon} .= $row[3]."\t".$row[4].";";
      $index{$chr}{$ind+1}{$id}{exon} .= $row[3]."\t".$row[4].";";
      $index{$chr}{$ind-1}{$id}{exon} .= $row[3]."\t".$row[4].";";
    }
  }
  close GFF;
  return %index;
}

sub exons {
  my ($self, $prefix, $genome) = @_;
  system ("gffread -O -w exons.fa -g ".$prefix."/reference/".$genome."_chr_all.fasta ".$prefix."/reference/".$genome."_genes.gff");
  open EXON, "<exons.fa" or die $!;
  my (%exon, $gene, $tran);
  while(my $exo = <EXON>){
  	chomp $exo;
  	if($exo =~ />(\w+)\.(\d+)/){
  		$gene = $1;
  		$tran = $2;
  	}elsif($exo=~ />(\w+)/){
  		$gene = $1;
  		$tran = 0;
  	}else{
  		$exon{$gene}{$tran} .= $exo;
  	}
  }
  close EXON;
  unlink "exons.fa";
  return %exon;
}

sub fas {
  my ($self, $prefix, $genome) = @_;
  my (%fas, $nam);
  open FAS, $prefix."/reference/".$genome."_chr_all.fasta" or die $!;
  while (my $line = <FAS>){
  	chomp $line;
  	if($line =~ />(.+)/){
      $nam = $1;
  		$fas{$nam} = "";
  	}else{
  		$fas{$nam} .= $line;
  	}
  }
  close FAS;
  return %fas;
}

sub geneinfo {
  my ($self, $prefix, $genome) = @_;
  my %geneinfo;
  open TRA, ">transcripts.fa" or die $!;
  my %exon = &exons($self, $prefix, $genome);
  foreach my $gene (sort keys %exon){
    my $long = 0;
    my ($longt, $longtr);
    foreach my $tran (keys %{$exon{$gene}}){
      if(length($exon{$gene}{$tran}) > $long){
        $long = length $exon{$gene}{$tran};
        $longt = $exon{$gene}{$tran};
        $longtr = $tran;
      }
    }
    $geneinfo{$gene} = $long;
    if($longtr){
      print TRA ">$gene\.$longtr\n$longt\n";
    }else{
      print TRA ">$gene\n$longt\n";
    }
  }
  close TRA;
  return %geneinfo;
}

sub mir {
  my ($self, $prefix, $genome) = @_;
  open GFF, "$prefix/reference/$genome\_miRNA_miRNA_star.gff" or die $!;
  my %mir;
  while (my $bb = <GFF>){
    chomp $bb;
    my @row = split /\t/, $bb;
    $mir{$row[8]}{"chromosome"} = $row[0];
    $mir{$row[8]}{"strand"} = $row[6];
    $mir{$row[8]}{"start"} = $row[3];
    $mir{$row[8]}{"end"} = $row[4];
  }
  close GFF;
  return %mir;
}

sub splitgff {
  my ($self, $prefix, $genome, $promoterLength) = @_;
  open GENE, "$prefix/reference/$genome"."_genes.gff" or die $!;
  open TE, "$prefix/reference/$genome"."_transposons.gff" or die $!;
  open TMP1, ">gene.gff" or die $!;
  open TMP2, ">promoter.gff" or die $!;
  while(my $aa = <GENE>){
    chomp $aa;
    my @row = split /\t/, $aa;
    if($row[2] =~ /gene/){
      if($row[8] =~ /ID=(\w+);/o){
        my $name = $1;
        if($row[8] =~ /Note=transposable_element_gene;/){
          $row[8] = $name."_TEG";
        }else{
          $row[8] = $name;
        }
        my $row = join "\t", @row;
        print TMP1 "$row\n";
        $row[8] .= "_promoter";
        if($row[6] eq "+"){
          $row[4] = $row[3] - 1;
          if($row[3] > $promoterLength){
            $row[3] = $row[3] - $promoterLength;
          }else{
            $row[3] = 1;
          }
        }else{
          $row[3] = $row[4] + 1;
          $row[4] = $row[4] + $promoterLength;
        }
        $row = join "\t", @row;
        print TMP2 "$row\n";
      }
    }
  }
  close GENE;
  close TMP1;
  close TMP2;
  open TMP, ">te.gff" or die $!;
  while(my $bb = <TE>){
    chomp $bb;
    my @row = split /\t/, $bb;
    if($row[2] =~ /transposable_element/){
      if($row[8] =~ /ID=(\w+);/o){
        $row[8] = $1;
        my $row = join "\t", @row;
        print TMP "$row\n";
      }
    }
  }
  close TMP;
  close TE;
}

sub ann {
  my ($self, $prefix, $genome, $binsize, $promoterLength) = @_;
  my %ann;
  if(-e $prefix."/reference/".$genome.".".$binsize.".annotation"){
  	open ANN, $prefix."/reference/".$genome.".".$binsize.".annotation" or die $!;
  		while (my $aa = <ANN>){
  			chomp $aa;
  			my @row = split /\t/, $aa;
  			my $id = shift @row;
  			$ann{$id} = join "\t", @row;
  		}
  	close ANN;
  }else{
    Ref->splitgff($prefix, $genome, $promoterLength);
    open GENE, "gene.gff" or die $!;
    open TE, "te.gff" or die $!;
    open PRO, "promoter.gff" or die $!;
    open MIR, $prefix."/reference/".$genome."_miRNA_miRNA_star.gff" or die $!;
    while (my $aa = <GENE>){
      chomp $aa;
      my @row = split /\t/, $aa;
      my $start = int($row[3] / $binsize);
      my $end = int($row[4] / $binsize);
      for(my $i=$start;$i<=$end;$i++){
        my $id = $row[0]."_".$i;
        $ann{$id} .= "GENE:".$row[8].";";
      }
    }
    while (my $bb = <TE>){
      chomp $bb;
      my @row = split /\t/, $bb;
      my $start = int($row[3] / $binsize);
      my $end = int($row[4] / $binsize);
      for(my $i=$start;$i<=$end;$i++){
        my $id = $row[0]."_".$i;
        $ann{$id} .= "TE:".$row[8].";";
      }
    }
    while (my $cc = <MIR>){
      chomp $cc;
      my @row = split /\t/, $cc;
      my $start = int($row[3] / $binsize);
      my $end = int($row[4] / $binsize);
      for(my $i=$start;$i<=$end;$i++){
        my $id = $row[0]."_".$i;
        $ann{$id} .= $row[8].";";
      }
    }
    while (my $dd = <PRO>){
      chomp $dd;
      my @row = split /\t/, $dd;
      my $start = int($row[3] / $binsize);
      my $end = int($row[4] / $binsize);
      for(my $i=$start;$i<=$end;$i++){
        my $id = $row[0]."_".$i;
        $ann{$id} .= "PROMOTER:".$row[8].";";
      }
    }
    close GENE;
    close TE;
    close MIR;
    close PRO;
    unlink ("gene.gff", "te.gff", "promoter.gff");

    foreach my $id (keys %ann){
      chop $ann{$id};
    }
  }
  return %ann;
}

sub lengthofchrom {
  my ($self, $prefix, $genome, $binsize) = @_;
  my %length;

  open FAI, $prefix."/reference/".$genome."_chr_all.fasta.fai" or die $!;
  while (my $fai = <FAI>){
    chomp $fai;
    my @row = split /\t/, $fai;
    $length{$row[0]} = int($row[1]/$binsize);
  }
  close FAI;
  return %length;
}

sub gann {
  my ($self, $prefix, $genome) = @_;
  my %gann = ();
  if(-e $prefix."/reference/".$genome.".BIN"){
    open BIN, $prefix."/reference/".$genome.".BIN";
    open FUN, $prefix."/reference/".$genome.".functional.annotation";
    while(my $aa = <BIN>){
      chomp $aa;
      my ($geneID, $bin) = split /\t/, $aa;
      $gann{$geneID} .= $bin.";";
    }
    close BIN;
    foreach my $gene (keys %gann){
      chop $gann{$gene};
    }
    while(my $bb = <FUN>){
      chomp $bb;
      my @fann = split /\t/, $bb;
      if($fann[0] =~ /(\w+)\.1$/){
        my $geneID = $1;
        if(exists $gann{$geneID}){
          $gann{$geneID} .= ",\"$fann[1]\",\"$fann[2]\",\"$fann[3]\",\"$fann[4]\"";
        }else{
          $gann{$geneID} .= "NA,\"$fann[1]\",\"$fann[2]\",\"$fann[3]\",\"$fann[4]\"";
        }
      }
    }
    close FUN;
  }
  return %gann;
}

sub PrimaryTranscript {
	my ($self, $prefix, $genome) = @_;
	system ("python ".$prefix."/scripts/getPrimaryTranscript.py ".$prefix."/reference/".$genome."_genes.gff > ".$genome.".PrimaryTranscript.txt");
	my (%primaryTranscript, $gene, $transcript, $exon);
	open PT, $genome.".PrimaryTranscript.txt" or die $!;
	while (my $aa = <PT>){
		chomp $aa;
		my @row = split /\t/, $aa;
		$primaryTranscript{$row[1]} = "";
	}
	close PT;
	open GTF, "<".$genome.".gtf" or die $!;
	open OUT, ">".$genome.".PrimaryTranscript.gtf" or die $!;
	while (my $bb = <GTF>){
		chomp $bb;
		my @row = split /\t/, $bb;
		if($row[8] =~ /transcript_id "(\w+\.\d+)"/){
			$transcript = $1; 
			print OUT "$bb\n" if(exists $primaryTranscript{$transcript});
		}
	}
	close GTF;
	close OUT;
	unlink $genome.".PrimaryTranscript.txt";
}

1;

__END__
