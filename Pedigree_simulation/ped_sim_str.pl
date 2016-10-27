
my $freq_file="STR_Freq_mutation";
my $locus_type="STR";

#$type is the pedigree scenario to test, could be SP (single parent), SS (single fullsib) or SH (single halfsib); 
#$group is case group, could be exp (relationship) or ctrl (unrelated control);
#$volume is the total number of cases to simulate, eg 100, 10000, ...
#my ($type,$group,$volume)=($ARGV[0],$ARGV[1],$ARGV[2]);
my ($type,$group,$volume)=("SP","exp",100);


my (%freq,%rate,%data);	#allele freq & mutation rate

open FREQ,"$freq_file";
while (<FREQ>) {
	if (!/\#/) {
		my @lines=split;
		$freq{$lines[0]}{$lines[1]}=$lines[2];
		$rate{$lines[0]}=$lines[3];
	}
}
close FREQ;

open OUT ,">result_$freq_file\_$type\_$group";
print OUT "\#round_id\tgroup\tperson1\tperson2\tlocus\tdata_person1\tdata_person2\tLR\n";

my @autoloci=keys %freq;

for my $round (1..$volume) {
	
	my (%grandpa,%grandma,%father,%mother,%uncle,%aunt,%cousin);
	my (@peoples,@h0_mothers,@h1_mothers,@h0_fathers,@h1_fathers,@known);
	my ($person1,$person2);
	
	if ($type eq "SP") {
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
	elsif ($type eq "SS") {
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
	
	elsif ($type eq "SH") {
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
		for my $i (1..$n) {	#initial the %pheno hash
			my $pheno_ith_people=$data{$peoples[$i-1]}{$locus};
			if ($pheno_ith_people eq "") {
				$pheno_ith_people="-,-";
			}
			my @temp=split /,/,$pheno_ith_people;
			my @allele1,@allele2;
			if ( $temp[0] eq "-" ) {
				@allele1=keys %{$freq{$locus}};
			}
			else {
				@allele1=($temp[0]);
			}
			if ( $temp[1] eq "-" ) {
				@allele2=keys %{$freq{$locus}};
			}
			else {
				@allele2=($temp[1]);
			}
			for my $ai (@allele1) {
				for my $aj (@allele2) {
					$pheno{$peoples[$i-1]}{$ai.",".$aj}++;
					if ($ai ne $aj) {
						$pheno{$peoples[$i-1]}{$aj.",".$ai}++;
					}
				}
			}
		}
		
		my @geno0=();
		my $prob_h0=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h0_fathers,\@h0_mothers,\@peoples,0 );
		my $prob_h1=&recurse( 0,\@geno0,\%pheno,$locus,$n,\@h1_fathers,\@h1_mothers,\@peoples,0 );
		my $subLR=$prob_h0/$prob_h1;
		print "$round\t$group\t$person1\t$person2\t$locus\t$data{$person1}{$locus}\t$data{$person2}{$locus}\t$subLR\n";
		print OUT "$round\t$group\t$person1\t$person2\t$locus\t$data{$person1}{$locus}\t$data{$person2}{$locus}\t$subLR\n";
	}
	
}


#####################################################################################

sub max {
	my @data=@_;
	my @temp=reverse sort {$a <=> $b} @data;
	return $temp[0];
}

sub min {
	my @data=@_;
	my @temp=sort {$a <=> $b} @data;
	return $temp[0];
}

sub randperson {
	my $male=$_[0];
	my %alleledata;
	
	for my $i (0..$#autoloci) {
		my @alleles=keys %{$freq{$autoloci[$i]}};
		my @genotypes;
		for my $j (0..1) {
			my $random=rand();
			for my $n (0..$#alleles) {
				if ($random <= $freq{$autoloci[$i]}{$alleles[$n]}) {
					$genotypes[$j]=$alleles[$n];
					last;
			}
				else {
					$random -= $freq{$autoloci[$i]}{$alleles[$n]};
				}
			}
		}
		$alleledata{$autoloci[$i]}=$genotypes[0].",".$genotypes[1];
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
		my @alleles=keys %{$freq{$autoloci[$i]}};
		my ($allele1,$allele2);
		my $rate1=$rate{$locus};
		my @p_types=split /,/,$father{$locus};
		my @m_types=split /,/,$mother{$locus};
		
		#paternal
		if (rand()<=0.5) {
			$allele1=$p_types[0];
		}
		else {
			$allele1=$p_types[1];
		}
		if (rand()<=$rate1) {
			if (rand()<=0.5) {
				$allele1+=1;
			}
			else {
				$allele1-=1;
			}
		}
		#maternal
		if (rand()<=0.5) {
			$allele2=$m_types[0];
		}
		else {
			$allele2=$m_types[1];
		}
		if (rand()<=$rate1/3.5) {
			if (rand()<=0.5) {
				$allele2+=1;
			}
			else {
				$allele2-=1;
			}
		}
		
		$child{$locus}=$allele1.",".$allele2;
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
	my $rate1=$rate{$locus};
	if ($k==$n) {
		return 1;
	}
	else {
		my %temp=%{$pheno{$peoples[$k]}};
		my @pheno_next=keys %temp;
		my $sum_p=0;
		my $f_index=$fathers[$k];
		my $m_index=$mothers[$k];
		my @geno_f=split /,/,$geno[$f_index-1];
		my @geno_m=split /,/,$geno[$m_index-1];
		
		my $f1=$geno_f[0];
		my $f2=$geno_f[1];
		my $m1=$geno_m[0];
		my $m2=$geno_m[1];
		my %tempf;
		my %tempm;
		$tempf{$f1}++;
		$tempf{$f2}++;
		$tempm{$m1}++;
		$tempm{$m2}++;
		
		for my $i (@pheno_next) {
			my ($tau_f,$tau_m)=(0,0);
			my @geno_i=split /,/,$i;
			my $mutation_step;
			if ($f_index>0) {
				$mutation_step=&min(abs($f1-$geno_i[0]),abs($f2-$geno_i[0]))+1E-4;
				if ($mutation_step<1E-3) {
					$tau_f=$tempf{$geno_i[0]}/2;
				}
				else {
					$tau_f=$rate1/(2*(10**( &max(1, int(abs($f1-$geno_i[0])+1E-4))-1 )))+$rate1/(2*(10**( &max(1, int(abs($f2-$geno_i[0])+1E-4))-1 )));
				}
			}
			else {
				$tau_f=&max(0.001,$freq{$locus}{$geno_i[0]});
			}
			if ($m_index>0) {
				$mutation_step=&min(abs($m1-$geno_i[1]),abs($m2-$geno_i[1]))+1E-4;
				if ($mutation_step<1E-3) {
					$tau_m=$tempm{$geno_i[1]}/2;
				}
				else {
					$tau_m=$rate1/(2*3.5*(10**( &max(1, int(abs($m1-$geno_i[1])+1E-4))-1 )))+$rate1/(2*3.5*(10**( &max(1, int(abs($m2-$geno_i[1])+1E-4))-1 )));
				}
			}
			else {
				$tau_m=&max(0.001,$freq{$locus}{$geno_i[1]});
			}
			my $subdel=$del+($tau_f<1E-12)+($tau_m<1E-12); 
			if ($tau_f<1E-20 && $tau_m<1E-20) {
				$sum_p+=0;
			}
			elsif ($subdel >=4) {	#To reduce computation load, high order small quantity is omitted.
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



	
