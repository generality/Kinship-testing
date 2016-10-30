#!/usr/bin/perl

my $data_file="demo_SNP_genotype";
my $freq_file="SNP_Freq_EAS";

my (%data,%freq,@loci,@alleles,%persons,%auto);
my $locus_type="SNP";

#$type is the pedigree scenario to test, could be SP (single parent), SS (single fullsib) or SH (single halfsib); 
#my $type=$ARGV[0];
my $type="SP";

my $mutation_rate=1.2E-8;
my $rare_allele_freq=1E-2;

if ($locus_type eq "SNP" ) {
	@alleles=("A","T","C","G");
}
open FREQ,"$freq_file";	#open freq file and store in %freq
while (<FREQ>) {
	my @lines=split;
	$freq{$lines[0]}{$lines[2]}=max($rare_allele_freq,$lines[3]);
	$freq{$lines[0]}{$lines[1]}=max($rare_allele_freq,1-$lines[3]);
}
close FREQ;

#open & get the PGM data (allele data of each locus)
open DATA,"$data_file";
my $titleline=<DATA>;
my @titles=split /\s+/,$titleline;
while (<DATA>) {
	my @lines=split /\s+/,$_;
	for my $i (1..$#lines) {
		if (&max(values %{$freq{ $titles[$i] }})>0) {
			$data{uc($lines[0])}{$titles[$i]}=$lines[$i];
			$persons{uc($lines[0])}=1;
			if (length($lines[$i])>1) {
				$auto{$titles[$i]}++;
			}
		}
	}
}
close DATA;
my @autoloci=keys %auto;

open OUT,">result_$freq_file\_$data_file\_$type";
print OUT "\#type\tperson1\tperson2\tlocus\tdata_person1\tdata_person2\tLR\n";

my @all_persons=sort {$a cmp $b}( keys %persons);

for my $p1 (0..$#all_persons) {
	my $person1=$all_persons[$p1];
	for my $p2 (0..$#all_persons) {
		my $person2=$all_persons[$p2];
		if ($p1 <= $p2) {
			next;
		}
		else {
			if ($type eq "SP") {
			#-----------------------------SP-------------------------------
				@peoples=("father","son");
				$peoples[1]=$person1;
				$peoples[0]=$person2;
				#h0
				@h0_mothers=(0,0);	#mother index of each people
				@h0_fathers=(0,1);	#father index of each people
				#h1
				@h1_mothers=(0,0);	#mother index of each people
				@h1_fathers=(0,0);	#father index of each people
			#--------------------------------------------------------------
			}
			
			elsif ($type eq "SS") {
			#-----------------------------SS-------------------------------
			@peoples=("father","mother","brother1","brother2");
				$peoples[3]=$person1;
				$peoples[2]=$person2;
				#h0
				@h0_mothers=(0,0,2,2);	#mother index of each people
				@h0_fathers=(0,0,1,1);	#father index of each people
				#h1
				@h1_mothers=(0,0,2,0);	#mother index of each people
				@h1_fathers=(0,0,1,0);	#father index of each people
			#--------------------------------------------------------------
			}
			
			elsif ($type eq "SH") {
				@peoples=("father","son1","son2");
				$peoples[2]=$person1;
				$peoples[1]=$person2;
				#h0
				@h0_mothers=(0,0,0);	#mother index of each people
				@h0_fathers=(0,1,1);	#father index of each people
				#h1
				@h1_mothers=(0,0,0);	#mother index of each people
				@h1_fathers=(0,1,0);	#father index of each people
			}
			
			#main function
			my $LR=1;
			my $n=$#peoples+1;	#people number
			for my $locus (@autoloci) {	#calculate LR for automosome loci
				my $subLR=1;
				my %pheno;
				
				for my $i (1..$n) {	#initial the genotype hash
					if ($locus_type eq "SNP" ) {
						my $pheno_ith_people=$data{$peoples[$i-1]}{$locus};
						if ($pheno_ith_people eq "") {
							$pheno_ith_people="NN";
						}
						my @temp=split //,$pheno_ith_people;
						my @allele1,@allele2;
						if ( ($temp[0] eq "N") || ($temp[0] eq "-") ) {
							@allele1=@alleles;
						}
						elsif (($temp[0] eq "A") || ($temp[0] eq "T") || ($temp[0] eq "C") || ($temp[0] eq "G") ) {
							@allele1=($temp[0]);
						}
						if ( ($temp[1] eq "N") || ($temp[1] eq "-") ) {
							@allele2=@alleles;
						}
						elsif (($temp[1] eq "A") || ($temp[1] eq "T") || ($temp[1] eq "C") || ($temp[1] eq "G") ) {
							@allele2=($temp[1]);
						}
						for my $ai (@allele1) {
							for my $aj (@allele2) {
								$pheno{$peoples[$i-1]}{$ai.$aj}++;
								if ($ai ne $aj) {
									$pheno{$peoples[$i-1]}{$aj.$ai}++;
								}
							}
						}
					}
				}
				my @geno0=();
				my $prob_h0=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h0_fathers,\@h0_mothers,\@peoples,0 );
				my $prob_h1=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h1_fathers,\@h1_mothers,\@peoples,0 );
				my $subLR=$prob_h0/$prob_h1;
				print OUT "$type\t$person1\t$person2\t$locus\t$data{$person1}{$locus}\t$data{$person2}{$locus}\t$subLR\n";
			}
		}
	}
}


###############################################################################################################
#                                                                                                             #
#                                               ¡ýsub function¡ý                                              #
#                                                                                                             #
###############################################################################################################
sub max {
	my @temp=@_;
	@temp=sort {$b<=>$a} @temp;
	return $temp[0];
}
sub recurse {
	my $k=$_[0];
	my @geno=@{$_[1]};
	my %pheno=%{$_[2]};
	my $locus=$_[3];
	my $n=$_[4];
	my @fathers=@{$_[5]};
	my @mothers=@{$_[6]};
	my @peoples=@{$_[7]};
	my $del=$_[8];
	if ($k==$n) {
		return 1;
	}
	else {
		my %temp=%{$pheno{$peoples[$k]}};
		my @pheno_next=keys %temp;
		my $sum_p=0;
		my $f_index=$fathers[$k];
		my $m_index=$mothers[$k];
		my $geno_f=$geno[$f_index-1];
		my $geno_m=$geno[$m_index-1];
		my $f1=substr($geno_f,0,1);
		my $f2=substr($geno_f,1,1);
		my $m1=substr($geno_m,0,1);
		my $m2=substr($geno_m,1,1);
		my %tempf;
		my %tempm;
		$tempf{$f1}++;
		$tempf{$f2}++;
		$tempm{$m1}++;
		$tempm{$m2}++;
		
		for my $i (@pheno_next) {
			my $tau_f,$tau_m;
			if ($f_index>0) {
				$tau_f=&max($mutation_rate,$tempf{substr($i,0,1)}/2);
			}
			else {
				$tau_f=max($rare_allele_freq,$freq{$locus}{substr($i,0,1)});
			}
			if ($m_index>0) {
				$tau_m=&max($mutation_rate,$tempm{substr($i,1,1)}/2);
			}
			else {
				$tau_m=max($rare_allele_freq,$freq{$locus}{substr($i,1,1)});
			}
			
			my $subdel=$del+($tau_f<0.0001)+($tau_m<0.00001);
			
			if ($tau_f<1E-2 && $tau_m<1E-2) {
				$sum_p+=0;
			}
			elsif ($subdel >=4) {
				$sum_p+=0;
			}
			else {
				my @subgeno=@geno;
				push @subgeno,$i;
				$sum_p+=($tau_f*$tau_m)*(&recurse($k+1,\@subgeno,\%pheno,$locus,$n,\@fathers,\@mothers,\@peoples,$subdel));
			}
		}
		return $sum_p;
	}
}
