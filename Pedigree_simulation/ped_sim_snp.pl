
my $freq_file="SNP_Freq_EAS";

my %data,%freq,@loci,@alleles;
my $locus_type="SNP";

#$kinship_type is the pedigree scenario to test, could be SP (single parent), SS (single fullsib) or SH (single halfsib). One could extend relationship by imitate the code; 
#$group is case group, could be exp (relationship) or ctrl (unrelated control);
#$volume is the total number of cases to simulate, eg 100, 10000, ...
#my ($kinship_type,$group,$volume)=($ARGV[0],$ARGV[1],$ARGV[2]);
my ($kinship_type,$group,$volume)=("SP","exp",100);

my $mutation_rate=1.2E-8;
my $rare_allele_freq=1E-2;

if ($locus_type eq "SNP" ) {
	@alleles=("A","T","C","G");
}
open FREQ,"$freq_file";	#open freq file and store in %freq
while (<FREQ>) {
	if (!/\#/) {
		my @lines=split;
		$freq{$lines[0]}{$lines[2]}=max($rare_allele_freq,$lines[3]);
		$freq{$lines[0]}{$lines[1]}=max($rare_allele_freq,1-$lines[3]);
	}
}
close FREQ;

my @autoloci=keys %freq;

open OUT,">result_$freq_file\_$kinship_type\_$group";
print OUT "\#round_id\tgroup\tperson1\tperson2\tlocus\tdata_person1\tdata_person2\tLR\n";

for my $round (1..$volume) {
	my %grandpa,%grandma,%father,%mother,%uncle,%aunt,%cousin;
	my @peoples,@h0_mothers,@h1_mothers,@h0_fathers,@h1_fathers;
	my $person1,$person2;
	
		if ($kinship_type eq "SP") {
		#-----------------------------SP-------------------------------
			@peoples=("father","son");
			$person1=$peoples[1];
			$person2=$peoples[0];
			@known=($person2);
			%father=randperson(1);
			if ($group eq "exp") {
				%mother=randperson(0);
				%son=trios(\%father,\%mother,1);
			}
			elsif ($group eq "ctrl") {
				%son=randperson(1);
			}
			$data{$peoples[0]}=\%father;
			$data{$peoples[1]}=\%son;
			#h0
			@h0_mothers=(0,0);	#mother index of each people
			@h0_fathers=(0,1);	#father index of each people
			#h1
			@h1_mothers=(0,0);	#mother index of each people
			@h1_fathers=(0,0);	#father index of each people
		#--------------------------------------------------------------
		}
		
		elsif ($kinship_type eq "SS") {
		#-----------------------------SS-------------------------------
		@peoples=("father","mother","brother1","brother2");
			$person1=$peoples[3];
			$person2=$peoples[2];
			@known=($person2);
			%grandma=randperson(0);
			%grandpa=randperson(1);
			%father=trios(\%grandpa,\%grandma,1);
			if ($group eq "exp") {
				%uncle=trios(\%grandpa,\%grandma,1);
			}
			elsif ($group eq "ctrl") {
				%uncle=randperson(1);
			}
			$data{$peoples[2]}=\%father;
			$data{$peoples[3]}=\%uncle;	
			#h0
			@h0_mothers=(0,0,2,2);	#mother index of each people
			@h0_fathers=(0,0,1,1);	#father index of each people
			#h1
			@h1_mothers=(0,0,2,0);	#mother index of each people
			@h1_fathers=(0,0,1,0);	#father index of each people
		#--------------------------------------------------------------
		}
		
		elsif ($kinship_type eq "SH") {
			@peoples=("father","son1","son2");
			$person1=$peoples[2];
			$person2=$peoples[1];
			@known=($person2);
			%father=randperson(1);
			%mother1=randperson(0);
			%mother2=randperson(0);
			%son1=trios(\%father,\%mother1,1);
			if ($group eq "exp") {
				%son2=trios(\%father,\%mother2,1);
			}
			elsif ($group eq "ctrl") {
				%son2=randperson(1);
			}
			$data{$peoples[1]}=\%son1;
			$data{$peoples[2]}=\%son2;	
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
		my $prob_h0=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h0_fathers,\@h0_mothers,\@peoples );
		my $prob_h1=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h1_fathers,\@h1_mothers,\@peoples );
		my $subLR=$prob_h0/$prob_h1;
		print "$round\t$group\t$person1\t$person2\t$locus\t$data{$person1}{$locus}\t$data{$person2}{$locus}\t$subLR\n";
		print OUT "$round\t$group\t$person1\t$person2\t$locus\t$data{$person1}{$locus}\t$data{$person2}{$locus}\t$subLR\n";
	}
}
close OUT;

#####################################################################################
sub max {
	my @data=@_;
	my @temp=reverse sort {$a <=> $b} @data;
	return $temp[0];
}
sub randperson {
	my $male=$_[0];
	my %alleledata;
	for my $i (0..$#autoloci) {
		if (1) {
			my @alleles;
			for my $j (0..1) {
				my $random=rand();
				if ($random<=$freq{$autoloci[$i]}{"A"}) {
					$alleles[$j]="A";
				}
				elsif ($random<=$freq{$autoloci[$i]}{"A"}+$freq{$autoloci[$i]}{"T"}) {
					$alleles[$j]="T";
				}
				elsif ($random<=$freq{$autoloci[$i]}{"A"}+$freq{$autoloci[$i]}{"T"}+$freq{$autoloci[$i]}{"C"}) {
					$alleles[$j]="C";
				}
				else {
					$alleles[$j]="G";
				}
			}
			$alleledata{$autoloci[$i]}=$alleles[0].$alleles[1];
		}
	}
	return %alleledata;
}
sub trios {
	my %father=%{$_[0]};
	my %mother=%{$_[1]};
	
	my $childsex=$_[2];
	my %child;
	my @loci=keys %father;
	for my $locus (@loci) {
		my ($allele1,$allele2);
		if (1) {
			my $random=rand();
			if ($random<=$mutation_rate) {
				if (rand()<=0.5) {
					if (rand()<=0.5) {
						$allele1="A";
					}
					else {
						$allele1="T";
					}
				}
				else {
					if (rand()<=0.5) {
						$allele1="C";
					}
					else {
						$allele1="G";
					}
				}
			}
			elsif ($random<=0.5) {
				$allele1=substr($father{$locus},0,1);
			}
			else {
				$allele1=substr($father{$locus},1,1);
			}
			my $random=rand();
			if ($random<=$mutation_rate) {
				if (rand()<=0.5) {
					if (rand()<=0.5) {
						$allele2="A";
					}
					else {
						$allele2="T";
					}
				}
				else {
					if (rand()<=0.5) {
						$allele2="C";
					}
					else {
						$allele2="G";
					}
				}
			}
			elsif ($random<=0.5) {
				$allele2=substr($mother{$locus},0,1);
			}
			else {
				$allele2=substr($mother{$locus},1,1);
			}
			$child{$locus}=$allele1.$allele2;
		}
		
	}
	return %child;
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
