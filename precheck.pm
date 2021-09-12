#!/usr/bin/env perl
package precheck;
use Modern::Perl;

my $prefix = $1 if($0 =~ /^(.+)\/.+$/);

sub dependencies {
  my ($self) = @_;

  print STDERR "\nChecking dependent software...\n";

  my $cutadapt = qx(cutadapt --version);
  my $samtools = qx(samtools --version);
  my $bowtie = qx(bowtie --version);
  my $bowtie2 = qx(bowtie2 --version);
  my $featureCounts = qx(featureCounts -v 2>&1);
  my $shortstack = qx(ShortStack -v);
  my $bedtools = qx(bedtools --version);
  my $R = qx(R --version);
  my $star = qx(STAR --version);
  my $clipper = qx(clipper -h 2>&1);
  my $fasterq_dump = qx(fasterq-dump -h);
  my $gffread = qx(gffread --version 2>&1);
  my $chip = qx(configureHomer.pl -list 2>&1|grep "+");
  my $sc = qx(umi_tools -v);
  my $deeptools = qx(deeptools --version 2>&1);
  my $bedgraphtobigwig = qx(bedGraphToBigWig 2>&1);

  if($cutadapt =~ /^\d\./){
  	print STDERR "cutadapt version $cutadapt";
  }else{
  	die "Please install cutadapt!\n";
  }
  if($samtools =~ /samtools (\d\.\d+)/){
  	if($1 lt "1"){
  		die "Please install samtools v1.x!\n";
  	}else{
  		print STDERR "samtools version $1\n";
  	}
  }else{
  	die "Please install samtools v1.x!\n";
  }
  if($bowtie =~ /bowtie.+ (\d\.\d\.\d)/){
  	print STDERR "bowtie version $1\n";
  }else{
  	die "Please install bowtie!\n";
  }
  if($bowtie2 =~ /bowtie2-(\d\.\d\.\d)/){
    print STDERR "bowtie2 version $1\n";
  }else{
  	die "Please install bowtie2!\n";
  }
  if($featureCounts =~ /v(\d\.\d+\.\d+)/){
    print STDERR "featureCounts version $1\n";
  }else{
    die "Please install featureCounts!\n";
  }
  if($shortstack =~ /version 3/){
  	print STDERR "$shortstack";
  }else{
  	die "Please install ShortStack!\n";
  }
  if($bedtools ne ""){
  	print STDERR "$bedtools";
  }else{
  	die "Please install bedtools!\n";
  }
  if($R =~ /R version (\d\.\d\.\d+)/){
  	print STDERR "R version $1\n";
  }else{
  	die "Please install R\n";
  }
  if($star =~ /(\d.+)/){
  	print STDERR "STAR version $1\n";
  }else{
  	die "Please install STAR!";
  }
  if($clipper =~ /^Usage/){
    print STDERR "CLIPper installed\n";
  }else{
    die "Please install CLIPper";
  }
  if($fasterq_dump =~ /fasterq-dump.+(\d\.\d+\.\d)/){
  	print STDERR "fasterq-dump version $1\n";
  }else{
  	die "Please install fasterq-dump!";
  }
  if($gffread =~ /(\d\.\d+\.\d)/){
    print STDERR "gffread version $1\n";
  }else{
    die "Please install gffread!";
  }
  if($deeptools =~ /(\d\.\d+\.\d)/){
    print STDERR "deeptools version $1\n";
  }else{
    die "Please install deeptools!";
  }
  if($chip =~ /tair10/){
    print STDERR "homer tair10 installed\n";
  }else{
    die "Please install the tair10 package of homer!";
  }
  if($sc =~ /(\d\.\d+\.\d)/){
    print STDERR "UMI-tools version $1\n";
  }else{
    die "Please install the UMI-tools!";
  }
  if($bedgraphtobigwig =~ /bedGraphToBigWig v (\d)/){
    print STDERR "bedGraphToBigWig version $1\n";
  }else{
    die "Please install UCSC-bedGraphToBigWig!";
  }
  system("Rscript --vanilla ".$prefix."/scripts/checkPackages.R");

  print STDERR "Precheck completed!\n\n";

}

1;
