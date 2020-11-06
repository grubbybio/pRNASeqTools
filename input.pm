#!/usr/bin/env perl
package input;

use Cwd 'abs_path';

sub run {

  my ($self, $input) = @_;

  my (@tags, @files, @par);
  my @inputtags = keys %$input;
  foreach $inputtag (@inputtags){
    die "Please do not name the group with a numeric initial" if ($inputtag =~ /^\d/);
    push @par, $inputtag;
    my @inputfiles = split /\+/, $$input{$inputtag};
    if($inputfiles[0] ne ""){
      if($inputfiles[0] !~ /^\d+$/){
        for(my $i=1;$i<=$#inputfiles+1;$i++){
          if($inputfiles[$i-1] !~ /^\//){
            if($inputfiles[$i-1] !~ /,/){
              if($inputfiles[$i-1] =~ /^~\/(.+)/){
                $inputfiles[$i-1] = $ENV{"HOME"}."/".$1;
              }else{
                $inputfiles[$i-1] = abs_path "../".$inputfiles[$i-1];
              }
            }else{
              my ($file1, $file2) = split /,/, $inputfiles[$i-1];
              if($file1 =~ /^~\/(.+)/){
                $file1 = $ENV{"HOME"}."/".$1;
              }else{
                $file1 = abs_path "../".$file1;
              }
              if($file2 =~ /^~\/(.+)/){
                $file2 = $ENV{"HOME"}."/".$1;
              }else{
                $file2 = abs_path "../".$file2;
              }
              $inputfiles[$i-1] = $file1.",".$file2;
            }
          }
          push @tags, $inputtag."_".$i;
        }
        push @files, @inputfiles;
        push @par, $#inputfiles+1;
      }else{
        push @par, $inputfiles[0];
        for(my $i=1;$i<=$inputfiles[0];$i++){
          push @tags, $inputtag."_".$i;
        }
      }
    }else{
      die "Please provide enough biological replicates\ntype \"-h\" for more information";
    }
  }
  return (\@tags, \@files, \@par);
}

1;

__END__
