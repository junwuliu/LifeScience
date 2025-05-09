#!/usr/bin/perl -w
use strict; use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($soft,$out);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$out,
	"i:s"=>\$soft,
) or &USAGE;
#&USAGE unless ($in and $out);
my %hash;
$/="\^";
my @title;
open IN,"zcat $soft|" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @content = split/\n/;
	my $samplename;
	if ($content[0] =~ /SAMPLE = (.*?$)/){
		$samplename = $1;
	}
	for (@content){
		if (/!Sample_title = (.*?$)/){
			$hash{$samplename}{'sampleTitle'} = $1;
			push @title,"sampleTitle";
		}
		if (/!Sample_characteristics_ch1 = (.*?): (.*?$)/){
			$hash{$samplename}{$1} = $2;
			push @title,$1;
		}
	}
	#last;
}
close IN;

open OUT,">",$out or die $!;
my %count;
my @uniq_title = grep { ++$count{ $_ } < 2; } @title;
@uniq_title = sort @uniq_title;
my $character = join("\t",@uniq_title);
print OUT "SampleName\t$character\n";
for my $sample(sort keys %hash){
	my @value;
	for my $char(@uniq_title){
		my $value = exists $hash{$sample}{$char} ? $hash{$sample}{$char} : "None";
		push @value,$value;
	}
	my $values = join("\t",@value);
	print OUT "$sample\t$values\n";
}
close OUT;


sub USAGE {
	my $usage=<<"USAGE";
	This script is used for extract meta infomation from GSE family soft
#--------
Program: $0
Contact:
#--------
USAGE
	print $usage;
	exit 1;
}
