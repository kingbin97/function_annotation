#!/usr/bin/perl

my $files_pfam_hmmer = shift @ARGV;
my $file_pfam_out = shift @ARGV;

open GENEID,$files_pfam_hmmer or die;
open OUT,">$file_pfam_out" or die;

my @col1_Gene_ID =();


my %gid = ();
my %goid = ();
my %goterm = ();
my %gocategory = ();

while (<GENEID>){
	chomp;
	my ($a,$b,$c,$d) = (split /\t/,$_)[0,1,2,3];
	# my $desc_string = join("\t",($b,$c,$d))
	push @col1_Gene_ID,$a if (not exists $goid{$a});
	push @{$goid{$a}},"[".$b."]";
	push @{$goterm{$a}},"[".$c."]";
	push @{$gocategory{$a}},"[".$d."]";
}

for my $a (@col1_Gene_ID){
	if(exists $goid{$a}){
		my $id_string = join(";",@{$goid{$a}});
		my $desc_string = join(";",@{$goterm{$a}});
		my $targ_string = join(";",@{$gocategory{$a}});
		print OUT "$a\t$id_string\t$desc_string\t$targ_string\n";
	}else{
		die;
	}
}

close GENEID;
close OUT;
