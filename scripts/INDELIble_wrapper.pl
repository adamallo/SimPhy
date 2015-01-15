#!/opt/local/bin/perl -w

# Quick perl script to simulate sequence evolution using INDELible along gene trees simulated by SimPhy.
# It doesn't support input indelible [BRANCHES] [TREE] [EVOLVE] and just one [MODEL]
# Diego Mallo 2015

use strict;
#GSL
use Math::GSL::RNG qw/:all/;
use Math::GSL::Randist qw/:all/;
use Math::GSL::CDF qw/:all/;
my $rng = Math::GSL::RNG->new();

#Config
my $w_dir;
my $input_file;
my $my_length;
my $length_command;
my $n_threads=0; #Maximum by default

my @files;
my $curr_sp=-1;
my $sp_counter=0;
my $sequence_counter=1;
my $filehandwrite;
my $filehandread;
my $filehandcontrol;
my $partitions;
my $evolves;
my $out_models;
my $temp_model;
my $model_id;
my $backup=$/;
my $locus;
my $id;
my $n_digits=0;
my $trees;
my $itree;

sub parse_sampling
{
	my $sampling_string="";
	my $distribution_code="";
	my $is_integer="";
	my $ret_value=0;
	
	($sampling_string)=@_;
	$sampling_string=~m/^(.*?):/;
	$distribution_code=uc($1);
	if ($sampling_string=~/:i$/i)
	{
		$is_integer=1;
	}
	else
	{
		$is_integer=0;
	}
	my @params=();
	while ($sampling_string=~m/([+-]?((\d+(\.\d*)?)|(\.\d+)))/g)
	{
		push(@params,$1);
	}
	if ($distribution_code eq "F")
	{
		scalar @params == 1 or die "Incorrect number of parameters for a fixed value. Sampling command \"$sampling_string\"";
		$ret_value=$params[0];
	}
	elsif ($distribution_code eq "N")
	{
		scalar @params == 2 or die "Incorrect number of parameters sampling a normal distribution. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_gaussian_ziggurat($rng->raw(),@params);
	}
	elsif ($distribution_code eq "L")
	{
		scalar @params == 2 or die "Incorrect number of parameters sampling a lognormal distribution. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_lognormal($rng->raw(),@params);
	}
	elsif ($distribution_code eq "U")
	{
		scalar @params == 2 or die "Incorrect number of parameters sampling an uniform distribution. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_flat($rng->raw(),@params);
	}
	elsif ($distribution_code eq "E")
	{
		scalar @params == 1 or die "Incorrect number of parameters sampling an exponential distribution. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_exponential($rng->raw(),@params);
	}
	elsif ($distribution_code eq "G")
	{
		scalar @params == 2 or die "Incorrect number of parameters sampling a gamma distribution. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_gamma($rng->raw(),@params);
	}
	elsif ($distribution_code eq "SL")
	{
		scalar @params == 3 or die "Incorrect number of parameters sampling a lognormal distribution multiplied by a constant. Sampling command \"$sampling_string\"";
		$ret_value=gsl_ran_lognormal($rng->raw(),$params[0],$params[1])*$params[2];
	}
	else
	{
		die "Unknown distribution code \"$distribution_code\" present in the sampling command \"$sampling_string\"";
	}
	
	if ($is_integer==1)
	{
		return int $ret_value;
	}
	else
	{
		return $ret_value;
	}
}

if ($#ARGV != 3)
{
	die "Incorrect number of parameters, Usage: script.pl directory input_config length numberofcores\n";
}

($w_dir,$input_file, $length_command,$n_threads)=@ARGV;
chdir($w_dir) or die "Error changing the working dir to $w_dir\n";
$w_dir=~m/([^\/]*).?$/;
$w_dir=$1;
opendir (my $dirs_handler, ".");

my @dirs = grep {-d "./$_" && ! /^\.{1,2}$/} readdir($dirs_handler);

#Working on the input file to be replicated
$/=undef;
open($filehandcontrol,$input_file) or die "Error opening the input config file $input_file";
my $content=<$filehandcontrol>;
close($filehandcontrol);
$/=$backup;

#Delete commented out text
$content=~s/\/\/.*?(?=\n)//g;
$content=~s/\/\*.*?\*\///sg;
$content=~s/\n+/\n/sg;

#TYPE Parsing
$content=~m/(\[TYPE\].*?)((?=\[[A-Z])|$)/s;
my $type=$1;

#MODEL Parsing
my @models;
while ($content=~m/(\[MODEL\].*?)((?=\[[A-Z])|$)/sg)
{
	push(@models,$1);
}
scalar @models > 1 and die "Unsuported number of [MODEL] in the input file (only one can be used)";

#PARTITIONS Parsing
my @partitions;
while ($content=~m/(\[PARTITIONS\].*?)((?=\[[A-Z])|$)/sg)
{
	push(@partitions,$1);
}
scalar @partitions > 0 and die "Unsuported [PARTITIONS] in the input file";

#EVOLVES Parsing
my @evolves;
while ($content=~m/(\[EVOLVE\].*?)((?=\[[A-Z])|$)/sg)
{
	push(@evolves,$1);
}
scalar @evolves > 0 and die "Unsuported [EVOLVE] in the input file";

#BRANCHES Parsing
my @branches;
while ($content=~m/(\[BRANCHES\].*?)((?=\[[A-Z])|$)/sg)
{
	push(@branches,$1);
}
scalar @branches > 0 and die "Unsuported [BRANCHES] in the input file";

#MAIN LOOP
foreach my  $dir (@dirs)
{
	$sp_counter=int($dir);
	print "\n\n\nTreating gene trees from the replicate $sp_counter\n";
	
	#Inside newdir
	chdir($dir) or die "Error changing the working dir\n";
	
	#Gene tree copy and modification 

	#INDELIBLE
	print "\t\nGenerating the INDELIBLE control.txt file\n";
	open($filehandwrite,">"."control.txt") or die "Error opening the file\n";
	
	@files=<g_trees*.trees>;
	$id='';
	$out_models='';
	$trees='';
	$partitions='';
	$evolves='[EVOLVE] ';
	
	foreach my $file (@files)
	{
		open($filehandread,$file) or die "Error opening the file $file\n";
		$file=~m/g_trees(\d*)\.trees/;
		$n_digits=length($1);
		$locus=int($1);
		$/="";
		$itree=<$filehandread>;
		chomp($itree);
		close($filehandread);

		$id=sprintf("%.*d",$n_digits,$locus);
		$temp_model=$models[0];
		$temp_model=~s/(\[MODEL\])\s*(.*?)\s*\n/$1 $2$id\n/;
		$model_id="$2$id";
		$temp_model=~s/\$\((.*?)\)/parse_sampling($1)/ge;
		$out_models.="$temp_model\n";
		$trees.=sprintf("\[TREE\] T%.*d %s\n",$n_digits,$locus,$itree);
		$my_length=parse_sampling($length_command);
		$partitions.=sprintf("\[PARTITIONS\] T%.*d \[T%.*d %s %s\]\n",$n_digits,$locus,$n_digits,$locus,$model_id,$my_length);
		$evolves.=sprintf("T%.*d 1 %.*d\n",$n_digits,$locus,$n_digits,$locus);
		
		$sequence_counter+=1;
	}
	
	print $filehandwrite $type,"\n",$out_models,"\n",$trees,"\n",$partitions,"\n",$evolves;
	close($filehandwrite);
	
	print "\tFile created\n";
	
	chdir("..");
	
}
my $is_parallel=`command -v parallel`;
if ($is_parallel && $n_threads!=0)
{
	print "Parallel sequence simulation\n";
	system("ls -d [0-9][0-9][0-9] | parallel -P $n_threads 'cd {} && indelible'");
}
else
{
	print "Sequential sequence simulation\n";
	foreach my $dir (@dirs)
	{
		$sp_counter=int($dir);
		print "\nSimulating sequences for replicate $sp_counter\n";	
		#Inside newdir
		chdir($dir) or die "Error changing the working dir\n";
		system("indelible");
		chdir("..")
	}
}
exit;


