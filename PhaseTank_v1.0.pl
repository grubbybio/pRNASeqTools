#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $version_num = "1.0";

my $help_message = Full_help($version_num);

#die if no any input info:
unless($ARGV[0]){
	die "$help_message" ;
}

#######################################DECLARE VARIABLES##########################################

my ($genome,$transcript_file,$ncrna_file,$miRNA_file,$degradome_file,$result_dir,$mixed_file);

my ($help,$version,$run_time,$now_time,$lib,$lib_num);

my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst);

my ($pred_table,$alignment,$cluster,$phasiRNA,$TAS_GENE,$miRNA_target_TAS,$phasiRNA_target,$cascade_file,$target_genes,$run_log);

my (@libs,%chr_seq,%chr_len,%ref_info,%aborted_reads,%reads_beg,%reads_seq,%phased_ratio,%phased_num,%phased_score,%tasiRNA_target);

my $normalization = 20000000;

my $PHASED_RATIO = 0.3;
my $PHASED_NUMBER= 4;
my $PHASED_ABUN = 100;#the total phasiRNA reads abundance should more than $PHASED_ABUN.
my $PHASED_ISLAND = 105;
my $PHASED_DRIFT = 2;
my $phased_size = 21;
my $extended_maplen = 80;
my $most_map_time = 5;
my $min_reads_abun =1; # minimum reads abundance kept in the mixed reads libraries.

my $asrp = 0.129;
my $RSRP = 1;
my $percent = 0.05;
my $min_phasiRNA_abun = 50; #minimum phasiRNAs abundance kept for phasiRNA target prediction. 

my $CALL_RSRP;
my $target_p;
my $target_m;

($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime();

$year += 1900;
$mon  += 1;
$run_time = sprintf("%d.%02d.%02d_%02d.%02d",$year,$mon,$day,$hour,$min);

$result_dir ="OUTPUT_$run_time";

$mixed_file = "Mixed_reads";

############################### OPTIONS #######################################################
#get options:
GetOptions ('genome=s' => \$genome,
			'cdna=s' => \$transcript_file,
			'lib=s' => \$lib,
			'filter=s' => \$ncrna_file,
			'miR=s' => \$miRNA_file,
			'degradome=s' => \$degradome_file,
			'target=s' => \$target_genes,
			
			'trigger_miRNA!' => \$target_m,
			'phasiRNA_target!' => \$target_p,
			
			'ratio=f' => \$PHASED_RATIO,
			'number=i' => \$PHASED_NUMBER,
			'abun=i' => \$PHASED_ABUN,

			'READ_abun=i' => \$min_reads_abun,
			'phasiRNA_abun=i' => \$min_phasiRNA_abun,
			'drift=i' => \$PHASED_DRIFT,
			'size=i' => \$phased_size,
			'nor=i' => \$normalization,
			'island=i' => \$PHASED_ISLAND,
			'extendLEN=i' => \$extended_maplen,
			'max_hits=i' => \$most_map_time,

			'per=f' => \$percent,
			'rsrp=f' => \$RSRP,
			'CALL_RSRP!' => \$CALL_RSRP,

			'dir=s' => \$result_dir,
			
			'version' => \$version,
			'help' => \$help);

# If help, print and quit:
if ($help){
	die "$help_message\n";
}

# If version, print version and quit:
if ($version){
	die "PhaseTank version $version_num.\n";
}

unless ($transcript_file or $genome){
	die "No transcript file or genome data, input transcript by '--cdna' or input genome data by '--genome'.\n";
}

unless ($lib){
	die "No small RNA sequencing data, input by '--lib'.\n";
}

@libs = split /\,/,$lib;

$lib_num = @libs;

##warning if no corresponding files:
if ($target_m){

	unless ($miRNA_file){
		die "No miRNA_file input for searching triggered miRNAs\n";
	}
	unless ($degradome_file){
		die "No degradome_file input for searching triggered miRNAs\n";
	}
 
}

if ($target_p){
	
	unless ($target_genes){
		die "No target genes for phasiRNA target prediction, input the target file by '--target'.\n";
	}
	
	unless ($degradome_file){
		die "no degradome_file input for searching targets of phasiRNAs\n";
	}
	
}

if ($CALL_RSRP){
	if ($genome){
		die "PhaseTank do not use genome data to predict the RSRP value. Please remove your genome data\n";
	}
	unless ($transcript_file){
		die "PhaseTank predict the RSRP value only using the whole cDNA data, input it by '--cdna'\n";
	}
}
#################################### SET OUTPUT FILES ###########################################

#make directory for the output files, which referred as its running time:

system "mkdir $result_dir";

my $path="$result_dir" . "\/";

$pred_table =  $path . "Pred_tab_$run_time";

$alignment = $path . "Align_$run_time";

$phasiRNA = $path . "PhasiRNA_$run_time";

$run_log = $path . "Run_log_$run_time";

$cluster = $path ."Excised_cluster_$run_time";

if ($target_m){
	$TAS_GENE = "TAS_gene.fa";
}

if ($target_m and $target_p){
#if ($target_m){
	$cascade_file = $path ."Cascades_$run_time";
}

#make directory for bowtie index files:
my $bowtie_index = "Bowtie-index";

#remove the previous files if it would be used again:
if (-d $bowtie_index){

	system "rm $bowtie_index/*";
	
}else{

	system "mkdir $bowtie_index";
}

if (-d $mixed_file){
	system "rm $mixed_file";
}

################################### MAIN ########################################################

#Print the input files :

$now_time = Watch_time();

Print_report("\n############################## PhaseTank version$version_num ##############################\n\nRun time: $now_time\n\n");

Print_report("#Input files:\n\n");

if ($transcript_file){

	Print_report( "\tRerference data: $transcript_file\n");

}elsif($genome){

	Print_report( "\tRerference data: $genome\n");
}

Print_report("\tRead file:\n");

my $n=0;
foreach my $read(@libs){

	$n++;

	Print_report("\t\t\($n\): $read\n");
}

if ($ncrna_file){

	Print_report("\tncRNA file: $ncrna_file\n");
}

if ($miRNA_file){

	Print_report("\tmiRNA file: $miRNA_file\n");
}

if ($degradome_file){

	Print_report("\tDegradome file: $degradome_file\n\n");
}

#Check the input files:

$now_time = Watch_time();

Print_report("##Begin to check files at $now_time..\n");

Check_input_files();

Print_report("\n##Finshed checking files\n");

#Normalize the reads abundance and merge them:

$now_time = Watch_time();

if ($lib_num > 1){
	
	Print_report("\n\n#Begin to merge libraries at $now_time..\n\n");

}else{

	Print_report("\n\n#Begin to normalize library at $now_time..\n\n");
}

Print_report("#Begin to merge libraries at $now_time..\n\n");

Merge_libraries();

if ($ncrna_file){

	$now_time = Watch_time();
	Print_report("\n#Begin to exclude other ncRNA mappable reads at $now_time..\n\n");

	Exclude_ncRNAs() ;
}

if ($genome){

	$now_time = Watch_time();
	Print_report("\n#Begin to excise RNA_producing cluster from genome data at $now_time\n\n");
	
	Produce_rna_cluster();
	
	Print_report("\n\t#Finish cluster finding \n");	
}

$now_time = Watch_time();
Print_report("\n#Begin to map reads to reference seq and parse the mapped reads data at $now_time\n");

Map_reads_Fill_hash($transcript_file,"transcript.fa",\%reads_beg,\%reads_seq);

$now_time = Watch_time();
Print_report("\n#Begin to parse reference information at $now_time\n");

Parse_ref_files();

Print_report("\n##Finished parsing files!\n\n");

#Activated by option '--CALL_RSRP':
if ($CALL_RSRP){

	$now_time = Watch_time();
	Print_report("\n#Begin to compute ASRP and RSRP_cut_off at $now_time\n\n");
	
	Compute_RSRP();
}

#Print out the important values used:
Print_report("\n##Important Values in PhaseTank:\n\n");
Print_report("\tLibrary_number: $lib_num\n\tASRP: $asrp\n\tRSRP: $RSRP\n");
Print_report("\tphased ratio: $PHASED_RATIO\n\tphased num: $PHASED_NUMBER\n");
Print_report("\tphased size: $phased_size nt\n\n");

$now_time = Watch_time();
Print_report("#Begin to process transcripts into core program at $now_time..\n\n");

Core_process();

Print_report("#Ready to OUTPUT the files:\n\n");

Output_files();

#Activated by option '--target_m':
if ($target_m){
	
	$now_time = Watch_time();
	Print_report("#Begin to find the triggered miRNAs by CleaveLand4 at $now_time..\n\n");

	$miRNA_target_TAS = $path . "miRNA_target_$run_time";
	Call_Cleaveland_miRNA_target();
}

#Activated by option '--target_P':
if ($target_p){
	
	$now_time = Watch_time();
	Print_report("#Begin to find the phasiRNAs targeted transcripts by CleaveLand4 at $now_time, phasiRNAs abundance >= $min_phasiRNA_abun..\n\n");
	
	$phasiRNA_target = $path . "PhasiRNA_target_$run_time";
	Call_Cleaveland_phasiRNA_target();
}

#To find the internal correlations for regulating networks:
if ($target_m and $target_p){

	$now_time = Watch_time();
	Print_report("#Begin to find regulatory cascades at $now_time..\n");

	Finding_cascades();
}

#To remove the Mixed_reads file merged by the reads libraries:
system `rm $mixed_file`;

$now_time = Watch_time();
Print_report("##Finished all processes in PhaseTank at $now_time!!\n\n#You can check the output files in directory './$result_dir/'.\n\n");

################################# HELP INFO #######################################################

sub Full_help {

    my($ver) = @_;
    my $help_info = "\nPhaseTank.pl version $ver

########################## PhaseTank_v1.0.pl ##########################

Copyright (c) 2014 Qingli Guo

Author:
	Qingli Guo, Northwest A&F University, email:guoql.karen\@gmail.com

Version:
	version 1.0; July 30 2014

#######################################################################
##Running PhaseTank from the following command: 

Usage: \$perl PhaseTank.pl --genome <genome_file.fa> --lib <read_file_list> [options]

Or \$ perl PhaseTank.pl --cdna <cdna_file.fa> --lib <read_file_list> [options] 

The followings are the detailed descriptions of the arguments and options in the use of PhaseTank:

Arguments: 

--genome <string>. Supply PhaseTank with genome sequence in FASTA format as reference sequences. Or --cdna <string>. Also could supply PhaseTank with cdna sequence in FASTA format as reference sequences.

--lib <string>. Supply PhaseTank with a comma-separated list of file(s) containing reads in FASTA format.

Options:
--filter <string>. Supply PhaseTank with FASTA format of other ncRNA sequences. It can help PhaseTank to exclude the reads mapped to other ncRNAs (e.g. tRNA, rRNA, snoRNA). 

--miR <string>. Supply PhaseTank with a list of miRNAs in FASTA format for miRNA-directed PHAS gene cleavage detection. This option will be ignored without ‘—trigger_miRNA’.

--degradome <string>. Supply PhaseTank with a set of degradome sequencing reads in FASTA format for phasiRNA targets prediction. 

--target <string>. Supply PhaseTank with a FASTA format file containing the interested genes, among which to search the phasiRNA targets.

--trigger_miRNA. Tell PhaseTank to detect miRNA-directed TAS cleavage. It is inactive by default.

--phasiRNA_target. Tell PhaseTank to predict the phasiRNA targeting genes. It is inactive by default.

--ratio <float>. Set phased ratio cutoff value. The default is 0.3.

--number <int>. Set phased number cutoff value. The default is 4.

--abun <int>. The total abundance of phased reads in the phasiRNA cluster. Default is 100. Note, the default normalization level is per twenty millions (20,000,000, can be changed by ‘--nor <int>’), thus the default abundance value of 100 here is equal to setting 5 of RPM (reads per million).

--READ_abun <int>. The minimum reads abundance to keep for PhaseTank prediction. Default is 1, which means if one read abundance is less than 1, it will be abandoned.

--phasiRNA_abun <int>. Minimum read abundance of phasiRNAs for target prediction. This option will be ignored without ‘—phasiRNA_target’.

--drift <int>. Maximum phased drift. The default is 2.

--size <int>. Length of phased reads. The default is 21.

--nor <int>. Tell PhaseTank the normalization level for the libraries. Default is 20,000,000.

--island <int>. That is the maximum separation distance of two phasiRNAs in each cluster. The default is 84.

--extendLEN <int>. The length on each side of siRNA cluster (or phasiRNA cluster) will be excised from the reference sequence. The default is 80.

-- max_hits <int>. Tell PhaseTank the ‘-m’ cutoff while using Bowtie (‘-m’ represent the maximum mapped hits to the reference, if goes out the value, the reads will be filtered out). The default is 5 here. Note that with this parameter changed, the prediction results may fluctuate slightly in big and small dataset due to a few reads may be removed in the big dataset. 

--per <float>. Within 0.01-1.00. The top percentage of RSRP value of sequences was put to the later program. The default is 0.05 (5%).

--rsrp <float>. The RSRP value for PhaseTank to filter the sequences. Default is 1.

--CALL_RSRP. Tell PhaseTank to estimate RSRP cutoff from the given reads libraries, which is set from the top 5% (default, can be changed by ‘—per <float>’) of RSRP value of sequences for the later processes. It is inactive by default. You could active this module by ‘--CALL_RSRP’ when you analyze other organisms (should use whole cDNA as input references). Or you also can use the default value instead.

--dir <string>. Set the directory in which PhaseTank will write its output files. The default is 'OUTPUT_run_time/'.

--help. Print the help message and quit.

--version. Print PhaseTank version number and quit.



Type \'PhaseTank.pl --help\' for full list of options

";
    return $help_info;
}

##################################### SUB PROGRAMS ##########################################################

sub Check_input_files{

	if ($transcript_file){
		FASTA_examiner_cdna($transcript_file);
	}
	
	if ($genome){
		FASTA_examiner($genome);
	}
	
	if ($ncrna_file){
		FASTA_examiner($ncrna_file);
	}
	
	foreach my $file (@libs){
		FASTA_examiner($file);
	}
	
	if ($degradome_file){
		FASTA_examiner($degradome_file);
	}
	
	if ($miRNA_file){
		FASTA_examiner($miRNA_file);
	}
	
	if ($degradome_file){
		FASTA_examiner($degradome_file);
	}
	
	if ($target_genes){
		 FASTA_examiner($target_genes);
	}
	
}


sub FASTA_examiner{
	
	my ($file) = @_;
	my $flag="no";
	
	open FILE,"<$file" || die "FAIL: $file does not exsist!\n";
    
	Print_report("\n\tChecking $file...\t");
	
	while(my $line = <FILE>){
	
		if ($line=~/^>/){
			$flag="yes";
			if ($line =~/\s/){
				$line=~s/ /_/g;
			}
		}
	}
	
	unless ($flag eq "yes"){
	
		die "ERROR in:$file is not in FASTA format.\n";
	}
	
	close(FILE);
}

sub FASTA_examiner_cdna{
	
	my ($file) = @_;
	my $flag="no";
	
	open FILE,"<$file" || die "FAIL: $file does not exsist!\n";
    
	Print_report("\n\tChecking $file...\n");
	
	while(my $line = <FILE>){
	
		if ($line=~/^>/){
			$flag="yes";
			
			if ($line =~/\s/){
				$line=~s/ /_/g;
			}
		}
	}

	unless ($flag eq "yes"){
	
		die "ERROR in:$file is not in FASTA format.\n";
	}
	
	close(FILE);
}


sub Merge_libraries{
	my %hash;
	
	foreach my $file(@libs){
		
		open LIB,"<$file";
		
		my $count=0;
		my $abun=0;
		my %tmp_hash=();
		Print_report("\t",$file,"\t\t");
		
		while( my $line1= <LIB>){
			my $line2=<LIB>;
			chomp ($line1,$line2);
			
			if($line1 =~/^>/){
				
				$count++;
				
				my @tmps=split /_x/,$line1;
				$abun+=$tmps[1];
				
				$tmp_hash{$line2}=$tmps[1];
			}
		}
		
		Print_report("lines:$count\t\tabun:$abun\n");

		$count=0;
		
		foreach my $seq( sort {$tmp_hash{$b} <=> $tmp_hash{$a}} keys %tmp_hash){
			
			my $normarlized_abun=int($tmp_hash{$seq}/($abun * $lib_num)* $normalization);
			
			$count++;
		
			if (defined $hash{$seq}){
			
				$hash{$seq}+=$normarlized_abun;
				
			}else{
			
				$hash{$seq}=$normarlized_abun;
			}
		}
		close(LIB);
	}

	my $count=0;
	open MIX,">$mixed_file";

	foreach my $seq (sort {$hash{$b} <=> $hash{$a}} keys %hash){
		my $abun = $hash{$seq};

		if ($abun >$min_reads_abun){
			$count++;
			print MIX ">t$count","_x$abun\n",$seq,"\n";
		}
		
	}
	close(MIX);
}


sub Parse_ref_files {
 
	my ($seq,$beg,$end);
	
	my $hash_name = "none";
	
	my $flag;
	open REF,"<$transcript_file";

	while (my $line=<REF>){

		chomp $line;
		
        if ($line=~/^>(\S+)\s+(\d+)\s+(\d+)$/){
           
			# to fill the hash and initialize the value for $hash_name and $seq;
			if ($hash_name ne "none"){

				$hash_name =~s/>//;
				
				$ref_info{$hash_name}{'seq'} = $seq;
				$ref_info{$hash_name}{'beg'} = $beg;
				$ref_info{$hash_name}{'end'} = $end;

				#remove the previous value:
				$hash_name = ();
                $seq = ();                
				$beg = ();
				$end = ();

			}
				
			$hash_name = $1;
			$beg=$2;
			$end=$3;

		}elsif($line=~/^>(\S+)$/){
			if ($hash_name ne "none"){
				$hash_name=~s/>//;
				$beg = 1;
				$end = length($seq);
				$ref_info{$hash_name}{'seq'}=$seq;
				$ref_info{$hash_name}{'beg'}=$beg;
				$ref_info{$hash_name}{'end'}=$end;

				#remove the previous value:
				$hash_name = ();
				$seq = ();                
				$beg = ();
				$end = ();
				
			}
			$flag="alone";
			$hash_name = $1;
			
		}elsif ($line !~ /^>/){

				$seq.=$line;
			}
    }
	
	#fill %ref_info with the info of last transcript;

	$hash_name =~s/>//;
	$ref_info{$hash_name}{'seq'} = $seq;
	
	if ($flag eq 'alone'){
	
		$ref_info{$hash_name}{'beg'} = 1;
		$ref_info{$hash_name}{'end'} = length($seq);
		
	}else{
	
		$ref_info{$hash_name}{'beg'} = $beg;
		$ref_info{$hash_name}{'end'} = $end;
	}

	close(REF);
	
	Print_report("\t$transcript_file has been parsed!\n\n");
	
}

sub Produce_rna_cluster{
	my ($name,$line_num,$width,@loci,%chr_reads_beg,%chr_reads_seq);
	
	my $id="chr.fa";
	
	open CHR,"<$genome";
	
	open CLUSTER,">$cluster";
	
	Print_report("\tParsing $genome ..\n\n");

	while (my $line=<CHR>){
		chomp $line;
	
		if ($line =~/^>(\S+)$/){

			$name=$1;

			$line_num=0;

		}else{
			$width=length($line) unless ($width);

			$line_num++;

			$chr_seq{$name}{$line_num}=$line;
				
			$chr_len{$name}+=length($line);
		}
	}
	close(CHR);
	
	Map_reads_Fill_hash($genome,$id,\%chr_reads_beg,\%chr_reads_seq);
	
	foreach my $chr (keys %chr_len){

		@loci=();

		foreach my $key ( sort{$chr_reads_beg{$chr}{$a} <=> $chr_reads_beg{$chr}{$b}} keys %{ $chr_reads_beg{$chr}}){
			my $pos = $chr_reads_beg{$chr}{$key};

			push @loci,$pos;
		}

		Define_cluster($chr,$width,\@loci);

	}
	close(CLUSTER);
	$transcript_file = $cluster;
}


sub Define_cluster{
	my ($chr,$width,$loci_ref)=@_;
	
	my ($i,$lim,$cluster_beg,$cluster_end,$cluster_len,$cluster_num,$print_line);
	
	$lim = @$loci_ref;
	
	$cluster_num=0;
	
	my @index=(); 
	
	Print_report("\n\tSearching and extracting clusters from $chr ..\t");
	for ($i = 0; $i<$lim; $i++){

		if ($i == 0){

			push @index, $i;
			next;
			
		}else{
		
			my $distance=$$loci_ref[$i]-$$loci_ref[$i-1];
			
			if ($distance<=100){
				push @index,$i;
			
			}else{
				
				my $num=$#index;
				
				if ($num<=3){
					
					@index=();
					push @index,$i;
					
					next;
				}
				my $start_read_pos = $$loci_ref[$index[0]];
				my $last_read_pos = $$loci_ref[$index[$num]];
				
				next if ($last_read_pos-$start_read_pos <=3*$phased_size);
				
				if ($start_read_pos>$extended_maplen){
					
					$cluster_beg = $start_read_pos-$extended_maplen;
				
				}else{
					
					$cluster_beg = 1;
				}

				if ($last_read_pos+$extended_maplen+30<$chr_len{$chr}){
					
					$cluster_end = $last_read_pos+$extended_maplen+30;
				
				}else{
					
					$cluster_end = $chr_len{$chr};
				}
				
				$cluster_len = $cluster_end-$cluster_beg+1;
				
				$cluster_num++;

				Extract_cluster_seq($chr,$cluster_beg,$cluster_end,$cluster_len,$cluster_num,$width,);
				
				@index=();
				
				push @index,$i;
			}
		}
	}

	Print_report("\tfound $cluster_num clusters\n\n");
}


sub Extract_cluster_seq{

	my ($chr,$beg,$end,$len,$cluster_num,$width) = @_;
	
	my($return_line,$cluster_seq,$remainder1,$remainder2,$remain,$seq1,$seq2,$seq3,$part1,$part3,$start_line,$end_line,$part2_beg,$part2_end,$step);
	
	$remainder1 = $beg % $width;

	$start_line = int($beg/$width) + 1;
	
	$seq1=$chr_seq{$chr}{$start_line};
	
	unless ($remainder1 == 0){
		
		$part1 = substr($seq1,$remainder1-1);
	
	}else{
		
		$part1=$seq1;
	}

	
	$remain=$len-length($part1);
	
	$step = int($remain/$width);
	
	$part2_beg=$start_line+1;
	$part2_end=$start_line+$step+1;
	
	for (my $line=$part2_beg;$line<$part2_end;$line++){
		
		$seq2.=$chr_seq{$chr}{$line};
	}
	
	unless ($seq2){
		
		$seq2="";
	}
	
	$remainder2 = $remain % $width;
	
	$seq3=$chr_seq{$chr}{$part2_end};
	
	unless($remainder2 ==0){
	
		$part3=substr($seq3,0,$remainder2);
		
	}else{
		$part3="";
	}

	$cluster_seq = $part1.$seq2.$part3;

	print CLUSTER ">$chr","_$cluster_num\t$beg\t$end\n$cluster_seq\n";

}


sub Map_reads_Fill_hash{

	my ($file,$id,$hash1_ref,$hash2_ref)=@_;
	
	Print_report("#Making bowtie index for $file..\n\n");
	
	my $bwt_ret1 = ` bowtie-build $file ./$bowtie_index/$id `;

	Print_report("#Mapping reads to $file using bowtie..\n\n");

	my $bwt_ret2 = ` bowtie -f -a -v 0 -m $most_map_time -x ./$bowtie_index/$id $mixed_file reads_mapped.bwt `;	

	my $map_count=0;

	open READ,"<reads_mapped.bwt" || die "FAIL:reads_mapped.bwt file does not exsist!\n";

	while (my $line = <READ>){

		chomp $line;
		my ($read_id,$str,$ref_id,$beg,$seq,$quality,$count) = split /\t/,$line ;

		$map_count++;
		
		my $key2 = "$read_id\t$str\t$map_count";
	
		$beg+=2 if ($str eq '-');	#plus 2 is to adjust mapped position from anti-sense strand to sense strand; 

		$$hash1_ref{$ref_id}{$key2} = $beg;
		$$hash2_ref{$ref_id}{$key2} = $seq;
	}
	close(READ);
	
	system `rm reads_mapped.bwt`;
}


sub Exclude_ncRNAs {
	
	my $ncrna_mapped="reads_mapped_ncrna.bwt";
	
	my $new_reads="New_mixed_reads.fa";
	
	Print_report("#Making bowtie index for $ncrna_file..\n\n");
	
	my $bwt_ret1 = ` bowtie-build $ncrna_file ./$bowtie_index/ncrna.fa `;

	Print_report("#Mapping reads to $ncrna_file using bowtie..\n\n");	

	my $bwt_ret2 = ` bowtie -f -a -v 0 ./$bowtie_index/ncrna.fa $mixed_file $ncrna_mapped`;
	
	open READ,"<$ncrna_mapped" || die "FAIL:reads_mapped_ncrna.bwt file does not exsist!\n";
	
	while (my $line =<READ>){

		chomp $line;
		
        my @tmps = split /\t/,$line;

		$aborted_reads{$tmps[0]}="yes";
	}
	
	open NEWREADS,">$new_reads";
	
	open MIXED,"<$mixed_file";
		
	while (my $line1 =<MIXED>){
		my $line2=<MIXED>;
		
		if ($line1=~/^>(\S+)$/){
		
			unless (defined $aborted_reads{$1}){
				
				print NEWREADS "$line1$line2"; 
			}
			
		}
	}

	close(MIXED);
	close(NEWREADS);
	
    system `rm $mixed_file`;
	
	$mixed_file=$new_reads;
	
	system `rm $ncrna_mapped`;
	
	Print_report("\nReads mapped to $ncrna_file has been excluded!\n\n");
}
	

sub Core_process {
	
	open ALIGN, ">$alignment" || die "\t\n$alignment can not be written in!\n";

	my $count=0;
	
	foreach my $ref_id(keys %ref_info){
		
		my @mapped_reads_plus =();
		my @mapped_reads_minus = ();
		my @mapped_reads = ();

		my %bins = ();

		my $abundance=0;        
		my $len=length ($ref_info{$ref_id}{'seq'});
		my $srp = 0 ;
		my $rsrp = 0;
	
		my $start = 0;
		my $end = $len-1;

#		we discard the transcripts (a: srp lower than the cutoff; b: mapped on ssRNA);		
		
		Search_mapped_reads ($ref_id,\$abundance,\@mapped_reads_plus,\@mapped_reads_minus,\@mapped_reads);

		$srp = $abundance / $len;

		next if ($srp == 0);

		$rsrp = sprintf("%.3f",log($srp/$asrp));

		$ref_info {$ref_id}{'rsrp'} = $rsrp;

		
		next if ($rsrp < $RSRP);

		#if only single strand have mappable reads, we discard it:
		my $mapped_plus = @mapped_reads_plus;
		
		my $mapped_minus = @mapped_reads_minus;
 
		next if ( $mapped_plus == 0 or $mapped_minus ==0);
		Assign_bins($start,$end,\%bins);
			
		Parse_phased_situation($ref_id,$len,$abundance,\$count,\@mapped_reads,\%bins);  

	}

	Print_report("\n\tPhaseTank has found $count PHAS genes!\n\n");
	close(ALIGN);
 
}


sub Search_mapped_reads {

	my ($ref_id,$abun_ref,$plus_str_ref,$minus_str_ref,$mapped_ref) = @_;
	
	foreach my $key( sort{$reads_beg{$ref_id}{$a} <=> $reads_beg{$ref_id}{$b} } keys %{ $reads_beg{$ref_id}}){
			
		push @$plus_str_ref, $key if ($key=~/\+/) ;
		push @$minus_str_ref, $key if ($key=~/\-/) ;
		
		push @$mapped_ref, $key;
			
		my @tmps=split /\t/,$key;
		
		$$abun_ref+= &Calculate_abun($tmps[0]);
	}

}


sub Assign_bins{

	my ($beg,$end,$bins_ref)=@_;
    
	for (my $i=$beg, my $bin_x=1; $i <= $end ; $i++){

		if ($bin_x <= $phased_size){
			$$bins_ref{$i}=$bin_x;
			$bin_x++;

		}else{			#when we count $phased_size bins, we count another $phased_size_bin;
			$bin_x=1;
			$$bins_ref{$i}=$bin_x;
			$bin_x++;

		}
	}
}


sub Parse_phased_situation{

	my ($ref_id,$len,$total_abun,$count_ref,$mapped_reads_ref,$bins_ref) = @_;
	
	my %bin_abun = ();
	my %bin_pos = ();
	
	my %bin_query_21nt = ();
    my %bin_pos_21nt = ();
	
	my %new_bin_pos_21nt;
	my %modified_bins;
	my %modified_bin_abun;
	my %modified_bin_abun_21nt;
	
	my @phased_pos = ();
	my @new_phased_pos = ();
	my @querys=();
	my @renew_phased_pos=();
		
	my $phased_ratio = 0;
	
	my $phased_abun = 0;

	my $modified_total_abun = 0;
	
    my $phased_num = 0;
	my $phased_score = 0;

	my ($most_abun,$more_abun,$max_bin,$more_bin);
	my ($read_len, $abun, $pos, $bin_x);
	my ($max_index,$flag);

	foreach my $key (@$mapped_reads_ref){
	
		$read_len = length ($reads_seq{$ref_id}{$key});
        my @tmps = split /\t/, $key;
		$abun = &Calculate_abun($tmps[0]);   
   		
		$pos = $reads_beg{$ref_id}{$key};
		$bin_x = $$bins_ref{$pos};	
		
		Filled_bin_info_hash($bin_x,$abun,$pos,$key,$read_len,\%bin_abun,\%bin_pos,\%bin_query_21nt,\%bin_pos_21nt);		
	}

	Find_most_abun_bin(\%bin_abun,\$most_abun,\$more_abun,\$max_bin,\$more_bin);
	
    if (defined $bin_pos_21nt{$max_bin}){
		
		if ($bin_pos_21nt {$max_bin} =~ /\,/){
			
			@phased_pos = split /\,/,$bin_pos_21nt{$max_bin};
			
			$phased_num = @phased_pos;
			
			if ($phased_num >= $PHASED_NUMBER){
				$flag = "once";
				
				$max_index = $#phased_pos;
				my $index_num;
				
				($phased_num,$index_num) = &Compute_phased_num($max_index,$PHASED_NUMBER,\@phased_pos);
				
				if ($index_num > 1){
					my @tmp_pos = ();
					
					foreach my $pos (@phased_pos){
						push @tmp_pos,$pos unless($pos=~/discarded/);
					}
					
					@new_phased_pos	= sort {$a <=> $b} @tmp_pos;				
					
					$max_index = $#new_phased_pos;
					
					&Compute_phased_num($max_index,$phased_num,\@new_phased_pos);
					
					$flag = "again"; 
				}
			}

		}else{
			$phased_num = 1;
		}
	}else{

		$phased_num = 0;
    }
	
    if ($phased_num >= $PHASED_NUMBER){
		
		if ($flag eq "again"){
			@phased_pos = ();
			@phased_pos = @new_phased_pos;
		}
		
		foreach my $pos (@phased_pos){
			push @renew_phased_pos,$pos unless($pos=~/discarded/);
		}
		
		my %poss=();
		
		foreach my $pos (sort {$a <=> $b} @renew_phased_pos){
				
			if (defined $poss{$pos}){
				next;
			}
			$poss{$pos}="done";
				
			my $query = $bin_query_21nt{$max_bin}{$pos};

			if ($query=~/\,/){
				
				my @tmps = split /\,/,$query;
				push @querys, $tmps[0];
				push @querys, $tmps[1];
			
			}else{
				push @querys, $query;
			}
		}
		
		my $max_index_new =	$#renew_phased_pos;
		
		my $start = $renew_phased_pos[0];
		my $end = $renew_phased_pos[$max_index_new];
		
		$ref_info{$ref_id}{'first_pos'} = $start;
		my $mRNA_beg = $ref_info{$ref_id}{'beg'};
		
		my $cluster_beg = $mRNA_beg + $start;
		my $cluster_end = $mRNA_beg + $end;
		
		$ref_info{$ref_id}{'cluster_region'} = $cluster_beg.":".$cluster_end;
		
		Assign_bins($start,$end,\%modified_bins);
		
		my $modified_total_abun = 0;
		
		foreach my $key (@$mapped_reads_ref){
			
			$pos = $reads_beg{$ref_id}{$key};
			
			$read_len = length ($reads_seq{$ref_id}{$key});
			
			if (defined $modified_bins{$pos}){
				
				$bin_x = $modified_bins{$pos};
				
				my @tmps = split /\t/, $key;
				
				$abun = &Calculate_abun($tmps[0]);
				$modified_total_abun += $abun;
				
				if (defined $modified_bin_abun{$bin_x}){

					$modified_bin_abun{$bin_x}+=$abun;		
				}else{
				
					$modified_bin_abun{$bin_x}=$abun;	
				}

				if  ($read_len == $phased_size){
					if (defined $modified_bin_abun_21nt{$bin_x}){
						$modified_bin_abun_21nt{$bin_x}+=$abun;		
					}else{

						$modified_bin_abun_21nt{$bin_x}=$abun;	
					}

					if (defined $new_bin_pos_21nt{$bin_x} ){

					$new_bin_pos_21nt{$bin_x}.=",$pos";
					
					}else{

						$new_bin_pos_21nt{$bin_x}=$pos;        
					}
				}
			}
		}
		
		Find_most_abun_bin(\%modified_bin_abun,\$most_abun,\$more_abun,\$max_bin,\$more_bin);
		
		if (defined $modified_bin_abun_21nt{$max_bin}){
			$phased_abun = $modified_bin_abun_21nt{$max_bin};
		}else{

			$phased_abun =0;
		}
		
		
		$phased_ratio = &Compute_phased_ratio($modified_total_abun,$most_abun,$more_abun,$max_bin,$more_bin);

		if ($phased_ratio >= $PHASED_RATIO and $phased_abun >= $PHASED_ABUN){
			
			$phased_ratio{$ref_id} = $phased_ratio;
			$phased_num{$ref_id} = $phased_num;
			
			$ref_info{$ref_id}{'phased_abun'} = $phased_abun;
			
			Output_phased_reads_map($ref_id,$len,$phased_abun,$phased_ratio,$phased_num,\@querys,\@renew_phased_pos);
			
			print ALIGN "\n##bin_x	Abundance	Abun_ratio(Abundance/total_abundance)	Phased_pos($phased_size nt mapped pos)\n";
			
			foreach my $key(sort {$modified_bin_abun{$b} <=> $modified_bin_abun{$a}} keys %modified_bin_abun){
				
				my $abun_ratio=sprintf ("%.3f",$modified_bin_abun{$key}/$modified_total_abun);

				if (defined $new_bin_pos_21nt{$key}){

					print ALIGN "bin_$key\t$modified_bin_abun{$key}\t$abun_ratio\t($new_bin_pos_21nt{$key})\n";

				}else{
				
					print ALIGN "bin_$key\t$modified_bin_abun{$key}\t$abun_ratio\t(_)\n";
				}
			}
			
			print ALIGN "\n\n";
			
			$phased_score = $phased_ratio * $phased_num * log ($phased_abun);

			$phased_score{$ref_id} = $phased_score;

			$$count_ref++;
			Print_report("\tpassed:\t$$count_ref\t$ref_id\n");
		}
	}
}


sub Calculate_abun{
    my ($read_id) = @_;
    my ($id,$abun)=split /_x/,$read_id;
    return $abun;
}


sub Filled_bin_info_hash{
	
	#we record the mapped pos into abun and pos hash one by one according to their bins;
	my ($bin_x,$abun,$pos,$key,$len,$bin_abun_ref,$bin_pos_ref,$bin_query_21nt_ref,$bin_pos_21nt_ref) = @_;
	my $count = 0;
	
	if (defined $$bin_abun_ref{$bin_x} and defined $$bin_pos_ref{$bin_x}){

		$$bin_abun_ref{$bin_x}+=$abun;
		$$bin_pos_ref{$bin_x}.= ",$pos";
		
	}else{

		$$bin_abun_ref{$bin_x}=$abun;
		$$bin_pos_ref{$bin_x}=$pos;
	}

	if  ($len == $phased_size){
		$count++;

		if (defined $$bin_pos_21nt_ref{$bin_x} ){

			$$bin_pos_21nt_ref{$bin_x}.=",$pos";
		
		}else{
           
		   $$bin_pos_21nt_ref{$bin_x}=$pos;        
		}
		
		if (defined $$bin_query_21nt_ref{$bin_x}{$pos}){
			
			$$bin_query_21nt_ref{$bin_x}{$pos}.=",$key";
		
		}else {
			
			$$bin_query_21nt_ref{$bin_x}{$pos}=$key;
		}		
	}	
}


sub Find_most_abun_bin {
	
	my ($bin_abun_ref,$most_abun_ref,$more_abun_ref,$max_bin_ref,$more_bin_ref)=@_;
	my $flag = "max";
    
	foreach my $key ( sort {$$bin_abun_ref{$b} <=> $$bin_abun_ref{$a}} keys %$bin_abun_ref ) {
    
		if ($flag eq 'max'){
		    
			$$most_abun_ref = $$bin_abun_ref{$key};
			$$max_bin_ref = $key;
			
			$flag = "second max";
			next;
			
		}
		
		if ($flag eq 'second max'){
		
			$$more_abun_ref = $$bin_abun_ref{$key};
			$$more_bin_ref = $key;
			
			my $flag= "end";
			last;
		}       
	}    
}


sub Compute_phased_ratio{

	my ($total_abun,$most_abun,$more_abun,$max_bin,$more_bin) = @_;

	my $phased_ratio = 0;  
	my $phased_abun_all_length = 0;
		
	my $fuzzy = abs($max_bin-$more_bin);
	my $max_drift = $phased_size - $PHASED_DRIFT;
	
	if ( $fuzzy <= $PHASED_DRIFT or $fuzzy >= $max_drift){	#whether the maximum bin and the second maximum bin are closed within 2 positions;
		
		$phased_abun_all_length = $most_abun+$more_abun;
		$phased_ratio= $phased_abun_all_length/$total_abun;
		
		}else{
		
		$phased_abun_all_length = $most_abun;
		$phased_ratio = $phased_abun_all_length/$total_abun;
	}
	
    return ($phased_ratio);
}

sub Compute_phased_num {
    
	my ($max_index,$filter,$phased_pos_ref)=@_;
    
	my @phased_nums;
	my $phased_number;

	my $start = $$phased_pos_ref[0];
    my $start_index = 0;
	
	my $passed_num = 1;
	
	for (my $index = 1 ; $index <= $max_index ; $index++){
        
		my $former_index = $index-1;
		my $distance = $$phased_pos_ref[$index]-$$phased_pos_ref[$former_index];
		
		if ($distance <= $PHASED_ISLAND){
			$passed_num++;

			if ($index == $max_index){
				if ($passed_num < $filter){
					for (my $i = $start_index ; $i <= $index ; $i++){
						$$phased_pos_ref[$i] = "discarded";
					}
				}else{
					push @phased_nums, $passed_num;
				}
            }
			
		}else{
			
			if( $passed_num < $filter){
				
				for (my $i = $start_index ; $i < $index ; $i++){
					$$phased_pos_ref[$i] = "discarded";                  
				}
				
			}else {
				push  @phased_nums, $passed_num;
			}
			
			if ($index == $max_index){
				$$phased_pos_ref[$index] = "discarded";
			}

			$start_index=$index;
			$passed_num = 1;
		}
    }
	  
	my @sorted_phased_nums = sort {$b <=> $a} @phased_nums;
	my $index_num = @sorted_phased_nums;
	
	if ($index_num == 0){
		$phased_number = 1;
	}else{
		$phased_number = $sorted_phased_nums[0];
	}

    return ($phased_number,$index_num);
}


sub Output_phased_reads_map{
	
	my ($ref_id,$len,$abun,$phased_ratio,$phased_num,$query_ref,$phased_pos_ref) = @_;

	my $ref_beg = $ref_info{$ref_id}{'beg'};
	my $ref_end = $ref_info{$ref_id}{'end'};
	
	my ($extended_beg,$extended_end,$new_beg,$new_end);
    my $end=$len-1;

	my $mapped_beg = $$phased_pos_ref[0];
	my $mapped_num = @$phased_pos_ref;
	my $mapped_end = $$phased_pos_ref[$mapped_num - 1] + $phased_size + $extended_maplen;
    
    my $mapped_len = $mapped_end - $mapped_beg + 1;

	if ($mapped_beg >= $extended_maplen){
		$extended_beg = $mapped_beg - $extended_maplen;
	}else {
		$extended_beg = 1;
	}
	
	if ($end <= $mapped_end) {
		$extended_end = $end;
	}else{
		$extended_end = $mapped_end;
	}
	
	$new_beg = $extended_beg + $ref_beg;
	$new_end = $extended_end + $ref_beg;
	
	$ref_info{$ref_id}{'excise_beg'} = $new_beg;
	$ref_info{$ref_id}{'excise_end'} = $new_end;
	
	my $excise_len = $extended_end - $extended_beg +1;
	$ref_info{$ref_id}{'excise_len'} = $excise_len;

	$phased_ratio = sprintf("%.3f", $phased_ratio);
	
	print ALIGN "//\n>$ref_id\texcise_beg:$new_beg\texcise_end:$new_end\texcise_len:$excise_len\tPhased_ratio:$phased_ratio\tPhased_num:$phased_num\tphased_abun:$abun\n";
	
	$ref_info{$ref_id}{'reads_query'} = ();
	
	foreach my $key (@$query_ref){		
		
		my $read_beg = $reads_beg{$ref_id}{$key};
		my $read_seq = $reads_seq{$ref_id}{$key};
		my $read_end = $read_beg + length($read_seq);
				
        if ($key=~/\+/){
			
			my $dot_num_before = $read_beg - $extended_beg;
			my $dot_num_after = $extended_end - $read_end;
			
			my @tmps = split /\t/, $key;
			
			print ALIGN "." x $dot_num_before,$read_seq,"." x $dot_num_after,"($tmps[0]\t$tmps[1]\t$read_beg)\n";
			
			my $read_abun = Calculate_abun($tmps[0]);
			
			if ($read_abun >= $min_phasiRNA_abun){
				$ref_info{$ref_id}{'reads_query'} .= ">". $ref_id . "_" . $read_beg . "(+)\n$read_seq\n";	
			}			
		}	
	}
	
	my $sense_seq = $ref_info{$ref_id}{'seq'};
	
	my $picked_sense_seq = substr($sense_seq,$extended_beg,$excise_len-1);

	$ref_info{$ref_id}{'excise_seq'} = $picked_sense_seq;
	
	my $antisense_seq = $picked_sense_seq;
	
	$antisense_seq =~tr/ATGC/TACG/;
		
	print ALIGN $picked_sense_seq,"(excised_seq\t$ref_id\t$extended_beg\t$extended_end)\n",$antisense_seq,"(anti-seq)\n";
	
	foreach my $key (@$query_ref){
		
		my $read_beg = $reads_beg{$ref_id}{$key}-2;
		my $read_seq = $reads_seq{$ref_id}{$key};
		
		if ($key =~ /\-/){
			
			$read_seq =~tr/ATGC/TACG/;
			
			my $dot_num_before = $read_beg-$extended_beg;
			my $dot_num_after = $extended_end-$read_beg - $phased_size;
			
			my @tmps = split /\t/, $key;
			print ALIGN "." x $dot_num_before,$read_seq,"." x $dot_num_after,"($tmps[0]\t$tmps[1]\t$read_beg)\n";
			
			my $read_abun = Calculate_abun($tmps[0]);
			
			if ($read_abun >= $min_phasiRNA_abun){
			
				$read_seq=reverse($read_seq);
				$ref_info{$ref_id}{'reads_query'} .= ">". $ref_id . "_" . $read_beg . "(-)\n$read_seq\n";
			}
		}
	}
}


sub Output_files {

	open TABLE, ">$pred_table";	
	print TABLE "#NO\tID\tBeg:End\tLength\tPhased_Ratio\tPhased_Abundance\tPhased_Number\tPhased_Score\n";
	
	if ($target_m){
		open TAS, ">$TAS_GENE";
	}
	
	open PHAS, ">$phasiRNA";
	
	my $num=0;
	
	foreach my $ref_id (sort {$phased_score{$b} <=> $phased_score{$a}} keys %phased_score){

		my $extended_beg = $ref_info{$ref_id}{'excise_beg'};
		my $extended_end = $ref_info{$ref_id}{'excise_end'};

		my $len = $ref_info{$ref_id}{'excise_len'};
		my $phased_abun = $ref_info{$ref_id}{'phased_abun'};
		my $ratio = $phased_ratio{$ref_id};

		my $phased_num = $phased_num{$ref_id};

		my $querys;
		
		if (defined $ref_info{$ref_id}{'reads_query'}){
			$querys = $ref_info{$ref_id}{'reads_query'};
		}
		
		my $seq = $ref_info {$ref_id}{'seq'};

		my $phase_score = $phased_score{$ref_id};
		
		$phase_score = sprintf("%.3f",$phase_score);
		
		$ratio = sprintf("%.3f", $ratio);
		
		$num++;

		my $table = "$num\t$ref_id\t$extended_beg\:$extended_end\t$len\t$ratio\t$phased_abun\t$phased_num\t$phase_score";
		
		$ref_info{$ref_id}{'table'} = $table;
		
		print TABLE "$table\n";
		
		if ($target_m){
		
			my $seq2 = $seq;
			$seq2 =~tr/ATGC/TACG/;
			
			$seq2=reverse($seq2);
			
			print TAS ">$ref_id\_1\n$seq\n>$ref_id\_2\n$seq2\n";
		}

		if (defined $ref_info{$ref_id}{'reads_query'}){
			print PHAS "$querys";
		}
	}

	close(TABLE);
	close(PHAS);

	if ($target_m){
		close(TAS);
	}

	Print_report("\t$alignment\n\t$pred_table\n\t$phasiRNA\n\n\thave been written completely!\n\n");
}


sub Call_Cleaveland_miRNA_target {
	my $siRNA='siRNAs.fa';
	system `cat $miRNA_file $phasiRNA > $siRNA`;
	system `perl ./CleaveLand4_modified.pl -e $degradome_file -u $siRNA -n $TAS_GENE -t -q > $miRNA_target_TAS`;

	Find_triggered_miRNA();

	open TAB,">$pred_table";
	
	my $line = "";
	my $count = 0;

	print TAB "#NO\tID\tBeg:End\tLength\tPhased_Ratio\tPhased_Abundance\tPhased_Number\tPhased_Score\tTriggered_miRNA\n";
	
	foreach my $ref_id (sort {$phased_score{$b} <=> $phased_score{$a}} keys %phased_score){
		
		if (defined $ref_info{$ref_id}{'triggered_miRNA'}){
			
			$line = $ref_info {$ref_id}{'table'} . "\t" . $ref_info {$ref_id}{'triggered_miRNA'};
			print TAB $line,"\n";
			$count++;
		
		}else{
			$line = $ref_info {$ref_id}{'table'} . "\t" . "NONE";
			print TAB $line,"\n";
		}
	}
	Print_report("\tThere are $count candidate TAS genes cleaved by input miRNAs\n\n");
	
	if ($line eq ""){
		Print_report("\tCan not find triggered miRNA for the candidate TAS gene!\n");
	}
	
	close(TAB);
	system `rm $TAS_GENE`;
	system `rm $siRNA`;
}


sub Call_Cleaveland_phasiRNA_target{
	
	$now_time = Watch_time();
	
	if ($target_genes){
		
		system `perl ./CleaveLand4_modified.pl -e $degradome_file -u $phasiRNA -n $target_genes -t -q > $phasiRNA_target`;
	
	}else{
		die "No target_genes input for phasiRNA target prediction, input by '--target',\n";
	}

}

	
sub Find_triggered_miRNA{
	
	open TARGET,"<$miRNA_target_TAS";

	while (my $line = <TARGET>){

		next if ($line =~ /^#/ or $line =~ /^SiteID/);
		
		chomp $line;
		my @tmps = split /\t/,$line;
		
		$tmps[2]=~/(\S+)\_(\d)$/;
		my $ref_id = $1;
		
		my $cleave_site = $tmps[5];
		
		my $phased_beg = $ref_info{$ref_id}{'first_pos'};
		if ($2 eq "2"){
		
			$cleave_site=length ($ref_info{$ref_id}{'seq'})-$cleave_site+2;
		}
		my $remainder = abs($cleave_site - $phased_beg) % $phased_size;

		if ($remainder >= 0 and $remainder <= 2 or $remainder==$phased_size-1 or $remainder == $phased_size-2){

			$ref_info{$ref_id}{'triggered_miRNA'}.= "$tmps[1]";
		}
	}
	close(TARGET);
}


sub Finding_cascades{
	
	open TARGET_M,"$miRNA_target_TAS";

	while (my $line =<TARGET_M>){
		
		if ($line =~/^#/ or $line =~/^Site/){
			next;
		}
		
		chomp $line;
		my @tmps = split /\t/,$line;
		
		unless (defined $ref_info{$tmps[2]}{'triggered_miRNA'}){
			$ref_info{$tmps[2]}{'cleaved'} = $tmps[1];
		}
	}
	close(TARGET_M);
	
	open TARGET_P, "$phasiRNA_target";
	
	while (my $line =<TARGET_P>){
		
		if ($line =~/^#/ or $line =~/^Site/){
			next;
		}
		
		chomp $line;
		
		my @tmps = split /\t/,$line;
		$tmps[1]=~/(\S+)\_\d+\(\S+\)/;
		my $name = $1;
		if (defined $ref_info{$name}{'phasiRNA'}){
			
			$ref_info{$name} {'phasiRNA'} .= "\t$tmps[1]";
		
		}else{
			$ref_info{$name} {'phasiRNA'} = $tmps[1];
		}
	
		if (defined $tasiRNA_target{$tmps[1]}){
		
			$tasiRNA_target{$tmps[1]} .= "\t$tmps[2]"; 
		
		}else{
			$tasiRNA_target{$tmps[1]} = $tmps[2];
		}
	
	}

	close(TARGET_P);
	
	open CASCADE ,">$cascade_file";
	
	foreach my $ref_id (sort {$phased_score{$b} <=> $phased_score{$a}} keys %phased_score){
		
		my $extended_beg = $ref_info{$ref_id}{'excise_beg'};
		my $extended_end = $ref_info{$ref_id}{'excise_end'};
		
		my $abun = $ref_info{$ref_id}{'phased_abun'};
		
		my $ratio = $phased_ratio{$ref_id};
		$ratio = sprintf("%.3f",$ratio);
		
		my $num = $phased_num{$ref_id};
	
		my $score = $phased_score{$ref_id};
	    $score = sprintf("%.3f",$score);
		
		print CASCADE "\n//\n>$ref_id\tExcise_region:$extended_beg:$extended_end\tPhased_ratio:$ratio\tPhased_abun:$abun\tPhased_num:$num\tPhased_score:$score\n";
		
		if (defined $ref_info{$ref_id}{'cleaved'}){
			
			my $cleaved_miRNA = $ref_info{$ref_id}{'cleaved'};
			
			print CASCADE "Cleaved by:\t$cleaved_miRNA\n";
		
		}elsif (defined $ref_info{$ref_id}{'triggered_miRNA'}){
			
			my $triggered_miRNA = $ref_info{$ref_id}{'triggered_miRNA'};
			
			print CASCADE "Triggered_miRNA:\t$triggered_miRNA\t-->\t$ref_id\n";	
		
		}

		if (defined $ref_info{$ref_id} {'phasiRNA'}){
			
			my %hash_exsit = ();
			
			print CASCADE "\nPhasiRNAs:\n";
			
			my @tmps = split /\t/,$ref_info{$ref_id} {'phasiRNA'};
			
			foreach my $phasiRNA (@tmps){
				
				if (defined $hash_exsit{$phasiRNA}){
					next;
					
				}else{
					
					$hash_exsit{$phasiRNA} = 'exist';
					print CASCADE "$phasiRNA\t-->\t$tasiRNA_target{$phasiRNA}\n";
				}
			
			}
		}

	}

	close(CASCADE);
	
}



sub Print_report{
	my ($info) = @_;
	
	print STDERR $info;
	
	open LOG,">>$run_log";
	print LOG $info;
	
}

sub Watch_time{

	($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime();
	
	$year += 1900;
	$mon  += 1;
	
	$now_time = sprintf("%d.%02d.%02d_%02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $now_time;
}


sub Compute_RSRP{
	
	my $count = 0;
	
	my $total_SRP = 0;
	
	my @srp_array = ();
#	$asrp = 0;
	
	foreach my $ref_id(keys %ref_info){
		
		my $abundance=0;        
		my $SRP = 0;

		$count++;
		Count_SRP($ref_id,\$SRP);

		$total_SRP += $SRP;
		push @srp_array, $SRP;
	}
	
	$asrp = $total_SRP / $count;
	$asrp = sprintf("%.3f",$asrp);	
	
	my @srp_new = sort {$b <=> $a} @srp_array;
	
	if ($percent > 0 and $percent <1){
		
		my $cut_number = int ($count * $percent);
		
		$RSRP = log($srp_new[$cut_number]/$asrp);
		
		$RSRP = sprintf ("%.3f",$RSRP);
		
		Print_report("count:$count\tpercent:$percent\tRSRP:$RSRP\n");
	
	}else{
		
		Print_report("--per [float 0..1], percent number should between 0 to 1.\n");
	}
	
}


sub Count_SRP {

	my ($ref_id,$SRP_ref) = @_;
	
	my $len=length ($ref_info{$ref_id}{'seq'});
	
	my $abun=0;
	
	foreach my $key( sort {$reads_beg{$ref_id}{$a} <=> $reads_beg{$ref_id}{$b} } keys %{ $reads_beg{$ref_id} }){
			
        my @tmps=split /\t/,$key;

        $abun += Calculate_abun($tmps[0]);
    }
	
	$$SRP_ref = $abun/$len;
}


##End of the sub-routines##

__END__

PhaseTank is a computational tool for genome-wide identification of phasiRNA involved regulatory cascades.

http://phasetank.sourceforge.net/

Cite: **.

----------------------------------------------------------------------------------------------
License

PhaseTank_v1.0.pl

Copyright (c) 2014 Qingli Guo

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

    PhaseTank: genome-wide computational identification of phasiRNAs and
    their regulatory cascades

AUTHOR
    Qingli Guo, Northwest A&F University, guoql.karen@gmail.com
	
   version 1.0; July 30 2014

-----------------------------------------------------------------------------------------------

PhaseTank can report detailed phased information of the predicted phasiRNA producing loci, and also the 

phasiRNAs regulatory cascades 'miRNA/phasiRNA -> PHAS gene -> phasiRNAs -> target genes'. 

To use PhaseTank, you need to install the related software (click 'Related Software' on the right of our website), and  in your PATH.
	
	1: perl (version 5.x)
	2: bowtie (version 0.12.x or 1.x)(Langmead, et al., 2009)
	
	The followings are required by cleaveland4 (Addo-Quaye, et al., 2009):
	3: Math::CDF (from CPAN)
	4: RNAplex (from Vienna RNA package)(Tafer and Hofacker, 2008)
	6: samtools(Li, et al., 2009)
	5: R

Besides, three perl scripts are also needed to add into your PATH:
	1: PhaseTank_v1.pl
	2: CleaveLand4_modified.pl
	3: GSTAr.pl

The PhaseTank package can be downloaded at 'Version 1.0', which includes a compress file named 'TUTORIAL_PhaseTank_V1.0.zip'. It includes:
	
	1: 1)	genome.fa (genome sequence, a multiline FASTA file; or cDNA sequence, a FASTA format file). It should be formatted as the list example.
	Argument : '--genome <string>'

example:	
>chr1
ATGTCCCTTCTGTTTCAACAGACAGTTCCTTTATCACACCTTCACAGGTCCCTCGATCCT
CCACTCTGCTTCCGCACTCACATACTGCTAATTCTTCTCCTGCTATCTCGACATCTTCCC
GGTTTCACAGGCTCTGATTGCGAATCTGCAGATCCTTCAATTGTCTCTGCGATTGCTCCT
GGAACTGCTACCACATCAGAAAGAGACTGTCCTGTGCGTACGGCAGGCTCAGATCCTGTT
CCTATTGGCGACAGCGGTACCTTTTTTGATGTTGGGACAGCTGCTCCTGAGCTACTTTCA
CCTAATAGACATCATATGATCACTCGGGCAAAGGATGGTATTCGCAAGCCTAATCCTCGT
TACAACCTGTTTACACAAAAATACACTCCCTCTGAACCAAAAACCATTACGTCTGCCTCC
CAGGATGGAGACAAGCTATGCAAGAAGAGATGTCGGCATTAA

	2: reads libraries 
	Argument : '--lib <file_1>,<file_2>,<file_n>' (comma-separated, NO SPACE between file_names is allowed )
	The small RNA-seq data, a multiline FASTA file. It should be formatted like:

example:
>t1_x4350713
TTTGGATTGAAGGGAGCTCTA
	
	(Note, 't1' is the id of the read, and can be any distinguished id name)
	
	The reason why we use mixed libraries is that the expression of TAS gene is dependent on the biological state of the tissue or cell.
	Thus, the merged data would provide us more information about this special class of RNA producing genes in a particular organism.
	While the user can also use one of their interested libraries to analyse the differential expression of PHAS genes in various genetic backgrounds.
		
	3: miRNAs (the microRNAs, a multiline FASTA file)
	Option : '--mi <string>'

example:	
>ath-miR156a
UGACAGAAGAGAGUGAGCAC
	
	It can be used to search the miRNA-triggered biogenesis for phasiRNAs.
	
	4: degradome data (a multiline FASTA file)
	Option : '--degradome <string>'

example:	
>read1
TTTTTTTTTTTTTTTTTTTT
>read2
TTTTTTTTTTTTTTTTTTTTT
>read3
TTTTTTTTTTTTTTTTTTTTT

	It can be used to validate the cleavage site for predicted targets by CleaveLand4.
	
	5: ncRNA file (the FASTA format seq of annotated other ncRNAs, a multiline FASTA file)
	Option : '--filter <string>'
	
example:

>ATMG01380.1
AAACCGGGCACTACGGTGAGACGTGAAAACACCCGATCCCATTCCGACCTCGATATGTGG
AATCGTCTTGCGCCATATGTACTGAGATTGTTCGGGAGACATGGTCCAAGCCCGGTGA

	6: target file (the FASTA format seq of candidate target for phasiRNA targets prediction)
	option: '--target <string>'
	
example:
>ATMG00010.1
ATGTCCCTTCTGTTTCAACAGACAGTTCCTTTATCACACCTTCACAGGTCCCTCGATCCT
CCACTCTGCTTCCGCACTCACATACTGCTAATTCTTCTCCTGCTATCTCGACATCTTCCC
GGTTTCACAGGCTCTGATTGCGAATCTGCAGATCCTTCAATTGTCTCTGCGATTGCTCCT
GGAACTGCTACCACATCAGAAAGAGACTGTCCTGTGCGTACGGCAGGCTCAGATCCTGTT
CCTATTGGCGACAGCGGTACCTTTTTTGATGTTGGGACAGCTGCTCCTGAGCTACTTTCA
CCTAATAGACATCATATGATCACTCGGGCAAAGGATGGTATTCGCAAGCCTAATCCTCGT
TACAACCTGTTTACACAAAAATACACTCCCTCTGAACCAAAAACCATTACGTCTGCCTCC
CAGGATGGAGACAAGCTATGCAAGAAGAGATGTCGGCATTAA
--------------------------------------------------------------------------------------------------------------------------

USAGE:
 
##Running PhaseTank from the following command: 

Usage: \$ perl PhaseTank.pl --genome <genome_file.fa> --lib <read_file_list> [options]
Or \$ perl PhaseTank.pl --cdna <cdna_file.fa> --lib <read_file_list> [options] 
The followings are the detailed descriptions of the arguments and options in the use of PhaseTank:

Arguments: 

--genome <string>. Supply PhaseTank with genome sequence in FASTA format as reference sequences. Or --cdna <string>. Also could supply PhaseTank with cdna sequence in FASTA format as reference sequences.

--lib <string>. Supply PhaseTank with a comma-separated list of file(s) containing reads in FASTA format.

Options:
--filter <string>. Supply PhaseTank with FASTA format of other ncRNA sequences. It can help PhaseTank to exclude the reads mapped to other ncRNAs (e.g. tRNA, rRNA, snoRNA). 

--miR <string>. Supply PhaseTank with a list of miRNAs in FASTA format for miRNA-directed PHAS gene cleavage detection. This option will be ignored without ‘—trigger_miRNA’.

--degradome <string>. Supply PhaseTank with a set of degradome sequencing reads in FASTA format for phasiRNA targets prediction. 

--target <string>. Supply PhaseTank with a FASTA format file containing the interested genes, among which to search the phasiRNA targets.

--trigger_miRNA. Tell PhaseTank to detect miRNA-directed TAS cleavage. It is inactive by default.

--phasiRNA_target. Tell PhaseTank to predict the phasiRNA targeting genes. It is inactive by default.

--ratio <float>. Set phased ratio cutoff value. The default is 0.3.

--number <int>. Set phased number cutoff value. The default is 4.

--abun <int>. The total abundance of phased reads in the phasiRNA cluster. Default is 100. Note, the default normalization level is per twenty millions (20,000,000, can be changed by ‘--nor <int>’), thus the default abundance value of 100 here is equal to setting 5 of RPM (reads per million).

--READ_abun <int>. The minimum reads abundance to keep for PhaseTank prediction. Default is 1, which means if one read abundance is less than 1, it will be abandoned.

--phasiRNA_abun <int>. Minimum read abundance of phasiRNAs for target prediction. This option will be ignored without ‘—phasiRNA_target’.

--drift <int>. Maximum phased drift. The default is 2.

--size <int>. Length of phased reads. The default is 21.

--nor <int>. Tell PhaseTank the normalization level for the libraries. Default is 20,000,000.

--island <int>. That is the maximum separation distance of two phasiRNAs in each cluster. The default is 84.

--extendLEN <int>. The length on each side of siRNA cluster (or phasiRNA cluster) will be excised from the reference sequence. The default is 80.

-- max_hits <int>. Tell PhaseTank the ‘-m’ cutoff while using Bowtie (‘-m’ represent the maximum mapped hits to the reference, if goes out the value, the reads will be filtered out). The default is 5 here. Note that with this parameter changed, the prediction results may fluctuate slightly in big and small dataset due to a few reads may be removed in the big dataset. 

--per <float>. Within 0.01-1.00. The top percentage of RSRP value of sequences was put to the later program. The default is 0.05 (5%).

--rsrp <float>. The RSRP value for PhaseTank to filter the sequences. Default is 1.

--CALL_RSRP. Tell PhaseTank to estimate RSRP cutoff from the given reads libraries, which is set from the top 5% (default, can be changed by ‘—per <float>’) of RSRP value of sequences for the later processes. It is inactive by default. You could active this module by ‘--CALL_RSRP’ when you analyze other organisms (should use whole cDNA as input references). Or you also can use the default value instead.

--dir <string>. Set the directory in which PhaseTank will write its output files. The default is 'OUTPUT_run_time/'.

--help. Print the help message and quit.

--version. Print PhaseTank version number and quit.

Type \'PhaseTank.pl --help\' for full list of options

----------------------------------------------------------------------------------------------------------

Analysis mode:
5.2 Analysis Demonstration for Different Modes
Here we provide four normal analysis modules using PhaseTank. In all of them, if you need to exclude some other ncRNAs you can input the FASTA format of the sequences by ‘--filter <ncRNA.fa>’ (for example: ‘--filter ath_ncRNA.fa’ in our datasets) to the following commands.

5.2.1	Predict PHAS loci from the given organism and read libraries

Irrelevant options: --miR, --degradome, --trigger_miRNA, --phasiRNA_target, --target, --phasiRNA_abun

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa 

5.2.2 Predict phasiRNAs and search their miRNA-triggered cleavage
Required options: --miR, --degradome, --trigger_miRNA
Irrelevant options: --target, --phasiRNA_target, --phasiRNA_abun

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --miR ath_miRNA.fa --degradome de_GSM278335.fa –trigger_miRNA

5.2.3 Predict phasiRNAs and their targets
Required options: --degradome, --target, --phasiRNA_target
Irrelevant options: --miR, --trigger_miRNA
Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --degradome de_GSM278335.fa --target ath_cDNA_TAIR10.fa --phasiRNA_target

5.2.4 Predict phasiRNAs, search the miRNA-triggered cleavage and detect phasiRNAs targets
	Required options: --miR, --degradome, --target, --trigger_miRNA, --phasiRNA_target

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --miR miRNA.fa --degradome de_GSM278335.fa --target ath_cDNA_TAIR10.fa --phasiRNA_target


Mode1: Predict phasiRNAs from the given organism and read libraries

Irrelevant options: --miR, --degradome, --trigger_miRNA, --phasiRNA_target, --target, --phasiRNA_abun

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa 

Mode2: Predict phasiRNAs and search their miRNA triggered cleaveage

Required options: --miR, --degradome, --trigger_miRNA
Irrelevant options: --target, --phasiRNA_target, --phasiRNA_abun

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --miR ath_miRNA.fa --degradome de_GSM278335.fa –trigger_miRNA

Mode3: Predict phasiRNAs and their targets

Required options: --degradome, --target, --phasiRNA_target
Irrelevant options: --miR, --trigger_miRNA
Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --degradome de_GSM278335.fa --target ath_cDNA_TAIR10.fa --phasiRNA_target

Mode4: Predict phasiRNAs, search the miRNA-triggered cleavage and detect phasiRNAs targets

Required options: --miR, --degradome, --target, --trigger_miRNA, --phasiRNA_target

Example:
\$ perl PhaseTank_v1.pl --genome ath_genome_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa --miR miRNA.fa --degradome de_GSM278335.fa --target ath_cDNA_TAIR10.fa --phasiRNA_target

------------------------------------------------------------------------------------------------------------------------------------

Output files will add to 'OUTPUT_runtime/'.	It includes:

	1: Pred_tab_runtime

	It is an tab-separated file and the "#" gives the what the value represent in each row.
	
	Note:
	
	If your analysis did not contain '--target_m', the last row 'Triggered_miRNA' will not exist.
	
	The 'No' is count by the phased_score. Bigger phased_score, smaller 'No' and more possible to be PHAS loci. 
	
	Beg:End(cluster) represent the begin and end position of phased cluster region, which is defined as a specific transcript region, which contained at least four phasiRNAs hits with a maximum separation distance of 84-nt based on the previous studies (Johnson, et al., 2009; Zhang, et al., 2012)
	
	2: Align_runtime
	
	Each PHAS loci info is separated by '//'.
	
	The '>' line is the loci information as marked.
	
	The following lines which contained many '.' are the mapped situation of each phased reads on each PHAS cluster. In the end, it contains the read mapped information.
	
	The cluster sequence is listed in the middle without any '.'. The top one is the sense strand and the bottom is anti-sense sequence.
	
	Reads over the cluster are mapped to the sense strand. In contrast, reads below the cluster are mapped to anti-sense strand.
	
	Lines begins with '##' are the annotation lines for the following information.
	
	The 'bin_x' is the bin number, which is from 'bin_1', 'bin_2' to 'bin_21' for each cluster see detail in the chapter written by Michael J. Axtell (2010). It ranked by their corresponding abundance.
	
	'Abun' is the total abundance of reads mapped to this bin position.
	
	'Abun_ratio' is the ratio of each bin position abun to the total abundance to this region.
	
	'Phased_pos(21 nt mapped pos)' are positions belonging to the bin_x and have mapped 21-nt reads. Like '1, 22, 43, 64, 85...155,176...260,281,302,...449' all these positions are belonging to 'bin_1'. But only part of the positions have 21nt mapped phased-reads.
	
	The same positions represent that phased reads are mapped on dsRNAs. We adjust the anti-sense mappable reads to the sense strand positions by plus 2 to the coordinates. Because the 21nt-duplex have 2-nt overhangs on the 3-terminal.
	
	3: miRNA_target_runtime
	
	This is the predicted results produced by CleaveLand4. The input file is the predicted PHAS loci in FASTA format, miRNAs in FASTA format and the degradome data.
	
	see detailed description in CleaveLand4_totorial.pdf (http://www.bio.psu.edu/people/faculty/Axtell/AxtellLab/Software.html)
	
	4: PhasiRNAs_runtime
	
	It is the FASTA format of phasiRNAs predicted by PhaseTank.

example:	
>AT2G27400.1_378(+)
TACAAGCGAATGAGTCATTCA

	Line initiates with '>' is the id of the phasiRNAs. For example, 'AT2G27400.1' is the ID for this phasiRNA producing gene. '378' is only a excision number to distinguish each phasiRNAs from the same loci. '(+)' is the strand it comes from.

	5: PhasiRNA_target_runtime
	
	This is the predicted results produced by CleaveLand4. The input file is the mRNA sequence in FASTA format, phasiRNAs in FASTA format and the degradome data.
	
	see detailed description in CleaveLand4_totorial.pdf (http://www.bio.psu.edu/people/faculty/Axtell/AxtellLab/Software.html)
	
	6: Cascades_runtime
	
	This is the cascades detected for predicted PHAS loci by PhaseTank in the given organism.
	
	Each PHAS loci is separated by '//'.
	
	Lines starts with '>' is the annotated information for this loci.
	
	The following line is:
	
	Sometimes, it reports triggered miRNAs for TAS genes, which is found to guide the cleavage of the TAS genes on the complementary region, especially on the 10-11 position  from the 5 terminal of the miRNA.
	
	Sometimes, it reports cleaved miRNAs for PHAS genes, which is the predicted results by CleaveLand4.
	
	If there is no targeting miRNAs, it directly reports the phasiRNAs and the predicted targets by CleaveLand4.
	
	If no targets by phasiRNA have been found, no information report.
	
	'AT2G39681.1_627(+)	-->	AT5G08680.1' represent the targeting direction from 'AT2G39681.1_627(+)' to 'AT5G08680.1'.
	
	7: Run_log_runtime
	it contains full list of the STEDRR screen print.
	
	8: Excised_cluster_runtime
	It contained all clusters excised from the genome.
example:
>chr_5_2	28430	28783
TTTGTTTTCGTGTTTTGTGGACTTAATTTGGGGGTTTATGATGAGTATGTGTAGGTATCTTTTTTTTTTTTTTTTTTTTTGTCAAAACACAACTTTCATTCATTAAGGCCTCAAGAGAGGAAAGTTGATACAAGCTACGATAATACAACAGAAAGAAAGATACAAGTTCCATGAGTTTTCCAAAGGGATGTTCACTTTTAAAAAGTTTGCAATGGTTGAAAAATCTGTTGAAACTATAGCGGATTCGACCATGAAGACGAGCTTCGCGAAATCCTCGACCAACGGGTGGTCATAGACGGCTCGAGCAAACAATGCATACAAGCTTTTAGAAACAAACTCCAAAGAGGAGGCTGT

"chr_5" is the chromosome number; "_2" is the number recorded in PhaseTank to distinguish each clusters excised from chr_5;
"28430" is the beg of this cluster; "28783" is the end of this cluster.
--------------------------------------------------------------------------------------------------------------

Conventions and recommendations

i.	In this manual, all the file names are in italic and the directory names are in bold and italic. 
	Besides, the command lines are listed in the grey backgrounds which start with $. Make sure your files are in UNIX format.

ii.	In our method, the relative small RNA production (RSRP) for a sequence. Therefore, the default value of RSRP here is set

	from the given data of Arabidopsis, which may be fluctuated in different datasets. However, the default value is recommended

	to use in your analysis. Because we have analyzed several libraries, the RSRP value will actually fluctuated for different 
 
	datasets, but it is quite slight. If you still want to estimate RSRP value in your dataset, you need to use whole cDNA sequences

	as your reference sequences and also should add option ‘--CALL_RSRP’ in your command line.

	For example:

	\$ perl PhaseTank_v1.pl --cdna ath_cdna_TAIR10.fa --lib GSM1174496.fa,GSM277608.fa,GSM342999.fa,GSM709567.fa,MTSRNA1.fa,RMMT10.fa –CALL_RSRP 

iii.	If there is a file containing other annotated ncRNAs in your aimed species, you can use this file to filter out the annotated ncRNAs in PhaseTank. 

iv.	In PhaseTank, the reference could be the genome sequences or any FASTA format of sequences (such as cDNA, EST, or your interested genes).

	According to the prediction results of our test, the genome sequences contained the richest information for PHAS genes detection.

	While if there is no complete genome assembly, other sequences data could also be used to predict PHAS loci with good sensitivity and specificity.


v.	The running time for PhaseTank mainly depends on the target prediction by CleaveLand4. It takes about 3-4 hours for analysis

	in 5.2.4 using the listed files and with the default settings. Option like ‘—phasiRNA_abun’ will largely influence the prediction time, because it directly decides the number of phasiRNAs which will be put into targets prediction pipelines.

vi.	We used modified CleaveLand4 in our pipeline for searching trigger miRNA and phasiRNA targets. The CleaveLand4 is just 

	modified to remove the screen output and some unimportant files, while the other parts and the core output file remain unchanged.

	Thus the prediction results will be clear.

--------------------------------------------------------------------------------------------------------------------------

This is the end of PhaseTank. Bugs, questions or suggestions are welcome for PhaseTank (Email: guoql.karen@gmail.com).
If you use PhaseTank in your work, please cite us :

PUBLISH
