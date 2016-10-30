my $datafile="demo.vcf";
my $idvfile="integrated_call_samples_v3.20130502.ALL.panel";

my $cutoff=0.03;

open DATA,"$datafile";
open OUT,">$datafile\_F_RT.txt";

open IN,"$idvfile";
my %idx;
<IN>;
while (<IN>) {
	my @lines=split;
	$idx{$lines[0]}=$lines[2];
}
close IN;

my $title=<DATA>;
my @titles=split /\s+/,$title;
my %cp;
$cp{"AFR"}=1;
$cp{"AMR"}=2;
$cp{"EAS"}=3;
$cp{"EUR"}=4;


print OUT "\#\#1\:\tAFR\n\#\#2\:\tAMR\n\#\#3\:\tEAS\n\#\#4\:\tEUR\n";
print OUT "\#CHROM\tPOS_ID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\_RT\n";

while (<DATA>) {
	my @lines=split;
	my @hcount;
	my @h1count;
	my @tcount;
	if ((length($lines[3])<2)&&(length($lines[4])<2)) {
		for my $i (9..$#lines) {
			my $popi=$cp{$idx{$titles[$i]}};
			if ($popi>0) {
				$tcount[$popi]++;
				if (($lines[$i]ne"0\|0")&&($lines[$i]ne"1\|1")) {
					$hcount[$popi]++;
				}
				elsif ($lines[$i]ne"0\|0") {
					$h1count[$popi]++;
				}
			}
		}
		$h1count[0]=&sum(@h1count);
		$tcount[0]=&sum(@tcount);
		$hcount[0]=&sum(@hcount);
		my @mafeach;
		my @heach;
		my $ft=($h1count[0]*2+$hcount[0])/($tcount[0]*2);
		$heach[0]=2*$ft*(1-$ft);
		$mafeach[0]=&min($ft,1-$ft);
		my $hs=0;
		for my $i (1..4) {
			my $ftemp=($h1count[$i]*2+$hcount[$i])/($tcount[$i]*2);
			$mafeach[$i]=&min($ftemp,1-$ftemp);
			$heach[$i]=2*$ftemp*(1-$ftemp);
			$hs+=($heach[$i]*$tcount[$i]);
		}
		$hs=$hs/$tcount[0];
		my $fst=1-$hs/$heach[0];
		if ($fst>$cutoff) {
			next;
		}
		for my $i (0..8) {
			print OUT "$lines[$i]\t";
		}
		print OUT "$fst\n";
	}
}


####################################################################
sub seqreverse {
	my $temp=$_[0];
	$temp=~s/A/W/gi;
	$temp=~s/T/A/gi;
	$temp=~s/W/T/gi;
	$temp=~s/C/W/gi;
	$temp=~s/G/C/gi;
	$temp=~s/W/G/gi;
	$temp =reverse $temp;				
}
sub mean {
	my @data=@_;
	my $n=0;
	my $sum=0;
	for (0..$#data) {
		$sum += $data[$_];
		$n++;
	}
	my $mean=$sum/$n;
}
sub min {
	my @data=@_;
	my @temp=sort {$a <=> $b} @data;
	$temp[0];
}
sub max {
	my @data=@_;
	my @temp=reverse sort {$a <=> $b} @data;
	$temp[0];
}
sub sum {
	my @data=@_;
	my $c = 0;
	foreach(@data){
		$c += $_;
	}
	return $c;
}
