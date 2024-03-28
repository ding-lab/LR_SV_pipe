#!/usr/bin/perl -w

#use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
 
&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This pipeline is designed for calling somatic variants from long-reads sequencing data.
        Author: ywzhang0713\@gmail.com Yuwei ZHANG 
        Usage: $0 config.txt

USAGE
print "$usage";
exit(1);
};


#print "Please be sure that the step1 LR_SV_pipe had been finished before using this script!";


my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;

my $thread = $conf{THREAD};
my $out = abs_path($conf{OUTDIR});
#make_path $out;
my $sample_f = abs_path($conf{SAMPLE});
#sample file should be in the below format:
#id sequencing_platform input_path

my $python = $conf{PYTHON};
my $pypy = $conf{PYPY};
my $samtools = $conf{SAMTOOLS};
my $seq_platform = $conf{PLATFORM};
my $ref = $conf{REF};
make_path abs_path($conf{OUTDIR});
my $bsub = "";
my $group = $conf{GROUP};
my $q = $conf{QUEUE};

my $cmd = "";
my $clair3 = "/storage1/fs1/dinglab/Active/Projects/yuweiz/software/ClairS/run_clairs";
my $clair_ot = "/storage1/fs1/dinglab/Active/Projects/yuweiz/software/ClairS-TO/run_clairs_to";

open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $input_number = scalar (@t);
	# print $input_number;
	my $id = shift @t;
	my $outdir = $out."/".$id."/Clair_output/";
	make_path $outdir;

	if($input_number == 3)
	{
		my $mode = "paired";
		my $t_input = abs_path($t[0]);
		my $n_input = abs_path($t[1]); 

		$cmd = "$clair3 --tumor_bam_fn $t_input --normal_bam_fn $n_input --ref_fn $ref --threads $thread  --output_dir $outdir --platform $seq_platform --output_dir $outdir --samtools $samtools --python $python --pypy $pypy --parallel /rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda3/envs/clair3/bin/parallel --clair3_path /rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda3/envs/clair3/bin/ --whatshap /rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda3/envs/clair3/bin/whatshap";
	}
	if($input_number == 2)
	{
		my $mode = "onlyTumor";
		my $t_input = abs_path($t[0]);
		$cmd = "$clair_ot --tumor_bam_fn $t_input --ref_fn $ref --threads $thread  --platform $seq_platform --output_dir $outdir --samtools $samtools --python $python --pypy $pypy --conda_prefix /rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda3/envs/clair3";
	}
	if(($input_number != 3) and ($input_number != 2))
	{
		$cmd="??";
		die("Please check your input sample list!");
	}


	open OUT, ">$out/$id/Clair_output/cmd.sh";
	print OUT "$cmd";
	close OUT;

	
	my $lsf_out = "$out/$id/bsub.log";
	my $lsf_err = "$out/$id/bsub.err";
	$bsub = "bsub -g $group -q $q -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 3000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $out/$id/Clair_output/cmd.sh";
	print $bsub."\n";
	system `$bsub`;
}
close IN;



sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}