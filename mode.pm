#!/usr/bin/env perl
package mode;
use Modern::Perl;
use MooseX::App;

opendir DIR, "." or die $!;
my $dir = join ",", readdir DIR;
closedir DIR;
has 'dir' => (
  is => 'ro',
  isa => 'Str',
  default => $dir,
  documentation => q[original files in the working directory.],
);

my $prefix = $1 if($0 =~ /^(.+)\/.+$/);
has 'prefix' => (
  is => 'ro',
  isa => 'Str',
  default => $prefix,
);

option 'outdir' => (
  is => 'rw',
  isa => 'Str',
  default => './out',
  documentation => q[Output directory.],
);
option 'genome' => (
  is => 'rw',
  isa => 'Str',
  default => 'ath',
  documentation => q[Currently supported genome: ath, osa, b73, gma, smo, bra, w22],
);
option 'thread' => (
  is => 'rw',
  isa => 'Int',
  default => '4',
  documentation => q[Threads used.],
);
option 'adaptor' => (
  is => 'rw',
  isa => 'Str',
  documentation => q[3' adaptor. 'AGATCGGAAGAGC' is the most common used, 'TGGAATTCTCGGG' is another.],
);
option 'control' => (
  is => 'rw',
  isa => 'HashRef',
  required => 1,
  documentation => q[seq files seperated by plus for control],
);
option 'treatment' => (
  is => 'rw',
  isa => 'HashRef',
  documentation => q[seq files seperated by plus for treatment, multiple treatments allowed],
);
option 'mask' => (
  is => 'rw',
  isa => 'Str',
  documentation => q[Masked sequences from the genome],
);

1;
__END__
