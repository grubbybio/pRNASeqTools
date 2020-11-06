#!/usr/bin/env perl
package validate_options;

use Modern::Perl;
use File::Path qw/make_path/;
use File::Copy qw/move/;

sub run {
  my ($self, $option) = @_;
  my %options=%$option;

  die 'Output directory exists! Please specify another output directory!' if(defined $options{'outdir'} and -e $options{'outdir'});
  make_path $options{'outdir'};
  move "log_".$main::time.".txt", $options{'outdir'}."/";
  chdir $options{'outdir'};
  if(defined $options{'adaptor'}){
    if($options{'adaptor'} eq "1"){
      $options{'adaptor'} = "AGATCGGAAGAGC";
    }elsif($options{'adaptor'} eq "2"){
      $options{'adaptor'} = "TGGAATTCTCGGG";
    }
  }
  if(defined $options{'adaptor2'}){
    if($options{'adaptor2'} eq "1"){
      $options{'adaptor2'} = "AGATCGGAAGAGC";
    }elsif($options{'adaptor2'} eq "2"){
      $options{'adaptor2'} = "TGGAATTCTCGGG";
    }
  }
  die 'Please use appropriate threads!' if(defined $options{'thread'} and $options{'thread'} =~ /[^1-9]/);
  die 'Please use an appropriate P value!' if(defined $options{'pvalue'} and $options{'pvalue'} > 0.05);
  die 'Please use an appropriate fold change!' if(defined $options{'foldchange'} and $options{'foldchange'} < 1.5);
  die 'Please use a supported strategy for mapping!' if(defined $options{'mmap'} and $options{'mmap'} =~ /[^ufrn]/);
  die 'Please use an appropriate window size!' if(defined $options{'binsize'} and $options{'binsize'} < 100);
  die 'Please specify an appropriate length of preferred small RNAs!' if(defined $options{'length'} and ($options{'length'} < 18 || $options{'length'} > 42));
  die 'Parameter conflict: nomapping and mappingonly!' if(defined $options{'nomapping'} and defined $options{'mappingonly'} and $options{'nomapping'} + $options{'mappingonly'} == 2);
  die 'Method not supported!' if(defined $options{'DESeq2Norm'} and $options{'DESeq2Norm'} ne 'DESeq2' and $options{'DESeq2Norm'} ne 'RPM');
  die 'Please specify an fasta file for mask!' if(defined $options{'mask'} and $options{'mask'} !~ /fasta$|fa$/);
  die 'Please select the correct style!' if(defined $options{'style'} and $options{'style'} !~ /histone|factor|tss/);
  die 'Cannot find the target list!' if(defined $options{'targets'} and $options{'targets'} ne "all" and !-e $options{'targets'});
  die 'Please use a suppported mode!' if(defined $options{'mode'} and $options{'mode'} ne 'mrna' and $options{'mode'} ne 'srna' and $options{'mode'} ne 'bulk' and $options{'mode'} ne 'sc');
  return %options;
}

1;

__END__
