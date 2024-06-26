#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw( min max sum);
use Scalar::Util qw(looks_like_number);

#Posiblle Databases
my $nr="NA";
my $interpro="NA";
my $kofam="NA";
my $nrDB="NA";
my $eggnog="NA";
my $deeploc="NA";
my $targetp="NA";
my $keggHierarchy="NA";

#Options
my $interproEvalue=0.00001;
my $nrEvalue=0.00001;
my $kofamEvalue=0.00001;
my $eggnogEvalue=0.00001;

#Out
my $out;

GetOptions(
	'-interpro=s' => \$interpro, #interproresults .tsv
	'-interproEvalue=s' => \$interproEvalue, #minimum evalue to keep interpro match [1e-5]
	'-nr=s' => \$nr, #bestmatch against nr .tsv
	'-nrDB=s'=> \$nrDB, #nr_24_08_2023/nr
	'-nrEvalue=s' => \$nrEvalue, #minimum evalue to keep nr match [1e-5]
	'-kofam=s' => \$kofam, #kofamScan results .tsv
	'-keggHierarchy=s' => \$keggHierarchy, #keggHierarchy database https://www.genome.jp/kegg-bin/get_htext : https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir= KO_Network.keg
	'-kofamEvalue=s' => \$kofamEvalue, #minimum evalue to keep Kofam match [1e-5]
	'-eggnog=s' => \$eggnog, #interproresults .tsv
	'-eggnogEvalue=s' => \$eggnogEvalue, #minimum evalue to keep eggnotch match [1e-5]
	'-deeploc=s' => \$deeploc, #deeploc .tsv
	'-targetp=s' => \$targetp, #targetp .tsv
	'-out=s' => \$out #prefix (.tsv will be added)
);

#Global variables
my $h; #Domain ID
my $desc; #domain description
my $score; # domain score
my @category; #list of database to print

#NR results
my %NRlist;
if ($nr ne "null"){
	push @category, "NR";
	print "NR results analysis...\n";
	open(NR,"$nr") or die "$nr not found";
	while (<NR>){
		chomp;
		my @cols=split("\t",$_);
		my ($gene,$frame)=$cols[0]=~/(.*)_(.)/;
		if(exists $cols[10] and $cols[10]<=$nrEvalue) {
			$h->{$gene}->{$frame}->{"NR"}="$cols[1];";
			$desc->{$gene}->{$frame}->{"NR"}="NA NA";
			$score->{$gene}->{$frame}->{"NR"}="$cols[10];";
			$NRlist{$cols[1]}="NA NA";
		}
	}
	close NR;
}

#Addition of NR taxonomy
if ($nrDB ne "NA" and $nr ne "NA"){
	print "NR database parsing...\n";
	open(NRDB,"$nrDB") or die "$nrDB not found";
	while (<NRDB>){
		chomp;
		if($_=~/>(.+?) (.+?)\[(.*?)\]/ and exists $NRlist{$1}){
			$NRlist{$1}="$3 $2";
		}
	}
	close NRDB;
}
	
#Kofamscan results
if ($kofam ne "NA"){
	push @category, "KO";
	print "KofamScan results analysis...\n";
	open(KO,"$kofam") or die "$kofam not found";
	while (<KO>){
		chomp;
		if ($_=~/^\*/){
			my $EC="";
			my($gene,$frame,$ko,$info,$evalue);
			if($_=~/\[(EC:.*)\]/){
				$EC=-$1;
				($gene,$frame,$ko,$evalue,$info)= $_=~ /([^ ]+?)_(.) +(K\d{5}) +[^ ]+ +[^ ]+ +([^ ]+) (.*) \[EC/ ;
			}
			else {
				($gene,$frame,$ko,$evalue,$info)= $_ =~ /([^ ]+?)_(.) +(K\d{5}) +[^ ]+ +[^ ]+ +([^ ]+) (.*)/ ;
			}
			unless (exists $h->{$gene}->{$frame}->{"KO"} and $h->{$gene}->{$frame}->{"KO"}=~/$ko$EC/){
				if($evalue<=$kofamEvalue){	
					$h->{$gene}->{$frame}->{"KO"}.="$ko$EC;";
					$desc->{$gene}->{$frame}->{"KO"}.="$info;";
					$score->{$gene}->{$frame}->{"KO"}.="$evalue;";
				}
			}
		}
	}
}
close KO;

#Kofam hierarchy
if ($keggHierarchy ne "NA" and $kofam ne "NA"){
	push @category, "BriteC";
	push @category, "BriteB";
	push @category, "BriteA";
	open(KH,"$keggHierarchy");

	my $kdes;
	my $kid;
	my %BriteId;
	my %BriteDesc;
	while(<KH>){
		chomp;	
		if($_=~/^(\D) *(K?\d+) +([^\[]+)/){
			my $rank=$1;
			my $id=$2;
			my $description=$3;
			if ($1 eq "D") {
				unless(exists $kid->{$id}->{"A"} and $kid->{$id}->{"A"}=~/$BriteId{"A"}/){
					$kid->{$id}->{"A"}.=$BriteId{"A"}.";";
					$kdes->{$id}->{"A"}.=$BriteDesc{"A"}.";";
				}
				unless(exists $kid->{$id}->{"B"} and $kid->{$id}->{"B"}=~/$BriteId{"B"}/){
					$kid->{$id}->{"B"}.=$BriteId{"B"}.";";
					$kdes->{$id}->{"B"}.=$BriteDesc{"B"}.";";
				}
				unless(exists $kid->{$id}->{"C"} and $kid->{$id}->{"C"}=~/$BriteId{"C"}/){
					$kid->{$id}->{"C"}.=$BriteId{"C"}.";";
					$kdes->{$id}->{"C"}.=$BriteDesc{"C"}.";";
				}
			}
			else{
				$BriteId{$rank}=$id;
				$BriteDesc{$rank}=$description;
			}
		}
	}
	close KH;
	#Addition of the Brite hierarchy in the table
	foreach my $gene (sort keys %$h){
		foreach my $frame ("P","1","2","3","4","5","6"){
			if (exists $h->{$gene}->{$frame}->{"KO"}){
				my (@kolist)=$h->{$gene}->{$frame}->{"KO"}=~/(K\d+)/g;
				foreach my $ko (@kolist){				
					foreach my $n ("A","B","C"){
						unless (exists $h->{$gene}->{$frame}->{"Brite$n"} and $h->{$gene}->{$frame}->{"Brite$n"}=~/$kid->{$ko}->{$n}/){
							$h->{$gene}->{$frame}->{"Brite$n"}.="$kid->{$ko}->{$n}";
							$desc->{$gene}->{$frame}->{"Brite$n"}.="$kdes->{$ko}->{$n}";
						}
					}
				}
			}		
		}
	}
}

#Interproscan results
if ($interpro ne "NA") {
	push @category, "IPR";
	push @category, "Pfam";

	open(IN,"$interpro") or die "$interpro not found";
	print "Interproscan results analysis...\n";
	while (<IN>){
		chomp;
		my @cols=split(/\t/,$_);
		my ($gene,$frame)=$cols[0]=~/(.*)_(.)/;
		unless(exists $h->{$gene}->{$frame}->{$cols[3]} and $h->{$gene}->{$frame}->{$cols[3]}=~/$cols[4]/){
			if($cols[8] eq "-" or $cols[8]<=$interproEvalue){
				$h->{$gene}->{$frame}->{$cols[3]}.="$cols[4];";
				$score->{$gene}->{$frame}->{$cols[3]}.="$cols[8];";
				$desc->{$gene}->{$frame}->{$cols[3]}.="$cols[5];";
			}
		}
		if (exists $cols[11] and $cols[11]=~/^IPR/){
			unless(exists $h->{$gene}->{$frame}->{"IPR"} and $h->{$gene}->{$frame}->{"IPR"}=~/$cols[11]/){
				$h->{$gene}->{$frame}->{"IPR"}.="$cols[11];";
				$desc->{$gene}->{$frame}->{"IPR"}.="$cols[12];";
				$score->{$gene}->{$frame}->{"IPR"}.="$cols[8];";
			}
		}
		if (exists $cols[13] and $cols[13]=~/^GO:/){
			my @GOs=split('\|',$cols[13]);
			foreach my $GO (@GOs) {
				unless(exists $h->{$gene}->{$frame}->{"GO"} and $h->{$gene}->{$frame}->{"GO"}=~/$GO/){
					$h->{$gene}->{$frame}->{"GO"}.="$GO;";
					$score->{$gene}->{$frame}->{"GO"}.="$cols[8];";
				}
			}
		}
	}
	close IN;

}

#eggNog results
if ($eggnog ne "NA"){
	push @category, "EGGNOG";
	print "EggNog results analysis...\n";
	open(EG,"$eggnog") or die "$eggnog not found";
	while (<EG>){
		if ($_=~/^#/){next;}
		chomp;
		my @cols=split(/\t/,$_);
		if ($cols[2]<=$eggnogEvalue){
			my ($gene,$frame)=$cols[0]=~/(.*)_(.)/;
			my ($id)=$cols[4]=~/(.+?)\@\d/;
			$h->{$gene}->{$frame}->{"EGGNOG"}="$id;";
			if($cols[8] eq "-"){$cols[8]="NA";}
			if($cols[9] eq "-"){$cols[9]="NA";}
			$desc->{$gene}->{$frame}->{"EGGNOG"}="$cols[8]|$cols[7]";
			$score->{$gene}->{$frame}->{"EGGNOG"}="$cols[2];";
		}
	}
}

#Deeploc
if ($deeploc ne "NA"){
	push @category, "DeepLoc";
	print "DeepLoc results analysis...\n";
	open(DL,"$deeploc") or die "$deeploc not found";
	my %header;
	while (<DL>){
		chomp;
		my @c=split(",",$_);
		if($c[0] eq "Protein_ID"){
			my $i=0;
			foreach my $k (@c){
				$header{$k}=$i;
				$i++;
			}
		}		
		else{
			my ($gene,$frame)=$c[0]=~/(.*)_(.)/;
			my @t=split("\\|",$c[1]);
			my $r;
			foreach my $v (@t) {
				$r.=1-sprintf("%.3f",$c[$header{$v}]).";"
			};
			$h->{$gene}->{$frame}->{"DeepLoc"}=$c[1];
			$desc->{$gene}->{$frame}->{"DeepLoc"}=$c[2];
			$score->{$gene}->{$frame}->{"DeepLoc"}="S=".substr($r,0,-1);
		}
	}
}

#TargetP
if ($targetp ne "NA"){
	push @category, "TargetP";
	print "TargetP results analysis...\n";
	open(TP,"$targetp") or die "$targetp not found";
	my %header;
	while (<TP>){
		chomp;
		my @c=split("\t",$_);
		if($c[0] eq "# ID"){
			my $i=0;
			foreach my $k (@c){
				$header{$k}=$i;
				$i++;
			}
		}		
		elsif ($c[0] !~ /^#/ and $c[1] ne "noTP") {
			my ($gene,$frame)=$c[0]=~/(.*)_(.)/;
			$h->{$gene}->{$frame}->{"TargetP"}=$c[1];
			$desc->{$gene}->{$frame}->{"TargetP"}=$c[7];
			$score->{$gene}->{$frame}->{"TargetP"}=(1-$c[$header{$c[1]}]);
		}
	}
}

if ($interpro ne "NA") {
	push @category, "TIGRFAM";
	push @category, "SMART";
	push @category, "Gene3D";
	push @category, "CDD";
	push @category, "GO";
}



#Results
print "Concatenation...\n";
open (OUT, ">$out.oneline.tsv");
open (OUT2, ">$out.long.tsv");
print OUT2 "Gene\tMethod\tFrame\tId\tDescription\tScore\n";
print OUT "Gene";

foreach my $cat (@category){
	print OUT "\t$cat frame|ID\t$cat desc\t$cat score";
}
print OUT "\n";
foreach my $gene (sort keys %$h){
	print OUT "$gene\t";
	foreach my $cat (@category){
		my $w="";
		my $x="";
		my $y="";
		my $z="1";
		foreach my $fr ("P","1","2","3","4","5","6"){
			#Keep the Protein predection except if a better score is available.
			my $NewScore;			
			if (exists $score->{$gene}->{$fr}->{$cat}){			
				$NewScore=$score->{$gene}->{$fr}->{$cat};
				if($NewScore=~/(.+?);/){
					$NewScore=$1;
				}
			}
			my $OldScore=$z;
			if($OldScore=~/(.+?);/){
				$OldScore=$1;
			}
			if ($fr eq "P" or (looks_like_number($z) and looks_like_number($NewScore) and $z > $NewScore)){
				if (exists $h->{$gene}->{$fr}->{$cat}) {
					($x) = $h->{$gene}->{$fr}->{$cat}=~/ *?(.*?);?$/;
				}
				if(exists $desc->{$gene}->{$fr}->{$cat}){
					($y) = $desc->{$gene}->{$fr}->{$cat}=~/ *?(.*?);?$/;
					if ($cat eq "NR" and $nrDB ne "NA" and exists $NRlist{$x}){
						$y = $NRlist{$x};
					}
					else {
						
					}
				}
				if(exists $score->{$gene}->{$fr}->{$cat}){
					($z) = $score->{$gene}->{$fr}->{$cat}=~/ *?(.*?);?$/;
				}
				$w=$fr;
			}
		}
		if($x ne ""){
			print OUT "$w|$x\t$y\t$z\t";
			print OUT2 "$gene\t$cat\t$w\t$x\t$y\t$z\n";
		}
		else { print OUT "\t\t\t";}
	}
	print OUT "\n";
}

close OUT;
close OUT2;
