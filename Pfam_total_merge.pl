#!/usr/bin/perl

my $files_pfam_hmmer = shift @ARGV;
my $file_pfam_out = shift @ARGV;

open GeneID,$files_pfam_hmmer or die;
open OUT,">$file_pfam_out" or die;

# my %h = ();
my @col1_Gene_ID =();
# my @col2_Pfam_ID =();
# my @col3_Pfam_Description = ();
# my @col4_Pfam_Target_Name = ();

my %gid = ();
my %pfamid = ();
my %pfamdesc = ();
my %pfamtarg = ();

while (<GeneID>){
	chomp;
	my ($a,$b,$c,$d) = (split /\t/,$_)[0,1,2,3];
	# my $desc_string = join("\t",($b,$c,$d))
	push @col1_Gene_ID,$a if (not exists $pfamid{$a});
	push @{$pfamid{$a}},"[".$b."]";
	push @{$pfamdesc{$a}},"[".$c."]";
	push @{$pfamtarg{$a}},"[".$d."]";
}

for my $a (@col1_Gene_ID){
	if(exists $pfamid{$a}){
		my $id_string = join(";",@{$pfamid{$a}});
		my $desc_string = join(";",@{$pfamdesc{$a}});
		my $targ_string = join(";",@{$pfamtarg{$a}});
		print OUT "$a\t$id_string\t$desc_string\t$targ_string\n";
	}else{
		die;
	}
}

close GeneID;
close OUT;
