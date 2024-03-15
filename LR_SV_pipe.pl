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

        This pipeline is designed for calling SVs from long-reads sequencing data.
        Author: ywzhang0713\@gmail.com Yuwei ZHANG 
        Usage: $0 config.tsv

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;

my $thread = $conf{THREAD};

my $out = abs_path($conf{OUTDIR});
my $sample_f = abs_path($conf{SAMPLE});
#sample file should be in the below format:
#id sequencing_platform input_path

my $step1 = "";
my $step2 = "";
my $step3 = "";
my $mode = $conf{SNIFFLES_MODE};
#### write you things ###
my $sample_size = 0;
my %allsample;


my $bsub = "";
my $group = $conf{GROUP};
my $q = $conf{QUEUE};

open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	my $pltform = $t[0];
	my $input = abs_path($t[1]);
	open OUT, ">$out/$id/cmd.sh";
	
	$step1 = "$conf{MINIMAP2} --MD -t $thread -ax $pltform $conf{INDEX} $input 2> $out/$id/00_mapping/$id.minimap2.log | $conf{SAMTOOLS} view -@ $thread -bS | $conf{SAMTOOLS} sort -@ $thread -o $out/$id/00_mapping/$id.minimap2.align.bam && $conf{SAMTOOLS} index -@ $thread $out/$id/00_mapping/$id.minimap2.align.bam";
	$step2 = "$conf{SNIFFLES} --input $out/$id/00_mapping/$id.minimap2.align.bam --vcf $out/$id/01_SVcall/$id.vcf.gz --snf $out/$id/01_SVcall/$id.snf 1> $out/$id/01_SVcall/$id.sniffles.log 2> $out/$id/01_SVcall/$id.sniffles.err";

	my $sniffles_version = `$conf{SNIFFLES} --version`;
	my $mosaic = "--mosaic";
	if(grep(/2.2/, $sniffles_version))
	{
		$mosaic = "--mosaic";
		# print "yes";
	}else{
		$mosaic = "--non-germline";
		# print "no";
	}
	$step3 = "$conf{SNIFFLES} --input $out/$id/00_mapping/$id.minimap2.align.bam --vcf $out/$id/02_SVcall_mosaic/$id.vcf.gz --snf $out/$id/02_SVcall_mosaic/$id.snf $mosaic 1> $out/$id/02_SVcall_mosaic/$id.sniffles.log 2> $out/$id/02_SVcall_mosaic/$id.sniffles.err";
	
	my $cmd = "";
	if($mode eq "Both")
	{
		make_path "$out/$id/00_mapping";
		make_path "$out/$id/01_SVcall";
		make_path "$out/$id/02_SVcall_mosaic";
		$cmd= $step1." && ".$step2." && ".$step3."\n";
	}
	if($mode eq "Basic")
	{
		make_path "$out/$id/00_mapping";
		make_path "$out/$id/01_SVcall";
		$cmd= $step1." && ".$step2."\n";
	}
	if($mode eq "Mosaic")
	{
		make_path "$out/$id/00_mapping";
		make_path "$out/$id/02_SVcall_mosaic";
		$cmd= $step1." && ".$step3."\n";
	}
	# print $cmd."\n";

	make_path abs_path($conf{OUTDIR});
	
	print OUT $cmd;
	close OUT;
	
	my $lsf_out = "$out/$id/bsub.log";
	my $lsf_err = "$out/$id/bsub.err";
	$bsub = "bsub -g $group -q $q -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $out/$id/cmd.sh";
	print $bsub."\n";
	system `$bsub`;

}
close IN;




# my $jobname = $conf{JOBNAME};
# my $partition = $conf{PARTITION};
# my $nodenum = $conf{NODENUM};
# my $task = $conf{NTASKS_PER_NODE};
# my $cpu = $conf{CPU_PER_TASK};
# my $cmdpath = "$out/cmd.sh";  #change
# my $failout = "$out/cmd.fail.sh"; #change
# my $array = $conf{ARRAY};

# open OUT, ">$out/slurm_cmd.sh";
# print OUT "\#\!\/bin\/bash
# \#SBATCH --job-name=$jobname
# \#SBATCH -p $partition
# \#SBATCH --nodes=$nodenum
# \#SBATCH --exclusive
# \#SBATCH --ntasks-per-node=$task
# \#SBATCH --cpus-per-task=$cpu
# \#SBATCH -o $out/$jobname.\%N.\%j.out 
# \#SBATCH -e $out/$jobname.\%N.\%j.err 
# \#SBATCH --array=$array

# arrayjob=\`cat $cmdpath | awk -v line=\$SLURM_ARRAY_TASK_ID \'{if (NR == line) print \$0}\'\`

# bash -c \"\$arrayjob \&\& {     echo \'Finish successfully\! \'; } || {     echo \$arrayjob >> $failout; }\"
# echo End time : \`date\`";

# close OUT;


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
