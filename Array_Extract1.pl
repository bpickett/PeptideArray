#!/usr/local/bin/perl -w
################################################################################################
##
## Perl Program to parse the results from .txt Peptide Array file, adjust for background, normalize, and summarize by median values
##
## Brett Pickett	July 2017
################################################################################################

use strict;
use warnings;

my $peptidePrefix = "JCVI ";#change this line for each array study
my $peptidePrefix1 = $peptidePrefix =~ s/\s//gr; #my $peptidePrefix1 = "JCVI";
my $peptideFilename = "JCVI-PrintedPeptideArrayList.tsv";

#Read all txt files from specified directory
my $PathToFolder = "."; #sets path to current directory
opendir (CURRENT, $PathToFolder);
my @filenames = readdir (CURRENT);
closedir (CURRENT);
#print "@filenames";

my @realFiles;
my ($basename, $suffix);
foreach my $filename1 (@filenames){ #iterates through all files in current directory
	#print $filename1;
	#Parse filename to get basename
	if ($filename1 =~m/(.+)[\.](.+)/){
		$basename = $1;
		$suffix = $2;
		if ($suffix eq 'txt'){
			#print $filename1;
			push (@realFiles, $filename1);
		}
	}
}

##read in mapping file (Sample_id -> Characterized virus)
my %characterizedSamples;
while (my $in = <DATA>){ #Reads text located after "__END__" in script
	chomp($in);
	if ($in =~m/(.*)?\t(.*)?/){
		my $sampleName = $1;
		my $virus = $2;
		$characterizedSamples{$sampleName} = $virus;
	}
}

##read in mapping file (peptide_id)
#format: >DENV1_E_120_134	KCVTKLEGKIVQYEN	 JCVI 1
#my %peptides;
my %name;
my %sequence;
my @sortedPeptides;
open (PEPTIDEFILE, "<$peptideFilename") || die "$peptideFilename: $!\n"; 	#check to see if file exists
	while (my $line1=<PEPTIDEFILE>) {
		chomp $line1;
		if ($line1=~m/^>(.*)\t(.*)\t$peptidePrefix(.*)/){
			my $pepLabel = $3;
#			$peptides{$pepLabel}{name} = $1;
#			$peptides{$pepLabel}{sequence} = $2;
			$name{$pepLabel} = $1;
			$sequence{$pepLabel} = $2;
			#print "label: $pepLabel\t name: $name{$pepLabel}\tsequence: $sequence{$pepLabel}\n";
			push @sortedPeptides, $pepLabel;
		}
	}
close (PEPTIDEFILE);
#print "@sortedPeptides";
my $sortedLength = @sortedPeptides;
#print "$sortedLength\n";

#my $testinput = "JCVI-3_H6_488H50P.txt";#delete
#push (@realFiles,$testinput);#delete
my %medianF_B;
#print "@realFiles";

#Process data from each file
my $sampleName;
print "Analyzing File: \n";
foreach my $realFilename (@realFiles){
	#print "$realFilename\n";
	if ($realFilename=~m/^$peptidePrefix1\-[0-9]+_([A-Z]+[0-9]+)_.*/){
		$sampleName = $1;
		print "\t...$sampleName\n";
	}
	open (ORIGFILE, "<$realFilename") || die "$realFilename: $!\n"; 	#check to see if file exists
	my $max_F_B = 0; #store maximum corrected median F-B value in array
	my %spotMedianF_B_storage;
	my $medianPeptideValue;
	while (my $line=<ORIGFILE>) {		#read all lines of the file
		my @results;
		if ($line=~m/^[0-9]+\t[0-9]+\t[0-9]+\t.*/){	#ignore the comment lines at the top of the csv file
			chomp $line;
			@results = split(/\t/,$line);
			if ($results[33] =~m/empty/ || $results[33] =~m/Flag[0-9]+/ || $results[33] =~m/HA[0-9]+/ || $results[33] eq "0") {#skip empty and control spots
				#print "skipping $results[33]\n";
				next;
				}
			else{
				my $spotMedianF = $results[8]; 
				my $spotMedianB = $results[14];
				my $spotName = $results[33];
				my $spotMedianF_B = $spotMedianF - $spotMedianB;
				#print "Name: $spotName: MedianF:$spotMedianF\tMedianB:$spotMedianB\tMedianF_B:$spotMedianF_B\n";
				
				#make all negative background-corrected values equal to zero
				if ($spotMedianF_B < 0){
					$spotMedianF_B = 0;
				}
				
				#identify brightest background-corrected non-control spot
				if ($spotMedianF_B > $max_F_B){
					$max_F_B = $spotMedianF_B;
					#print "NEW MAX VALUE: $max_F_B = $spotMedianF_B\n";
				}
				
				#capture peptide number and corrected value for each spot and store in hash
				if ($spotName =~m/$peptidePrefix([0-9]+)/){
					my $peptideNumber = $1;
					#my $a = 1;
					if (!exists $spotMedianF_B_storage{$peptideNumber}{1}){
						$spotMedianF_B_storage{$peptideNumber}{1} = $spotMedianF_B;
						#print "stored peptide #: $peptideNumber\t1\t$spotMedianF_B_storage{$peptideNumber}{1}\n";
					}
					elsif (!exists $spotMedianF_B_storage{$peptideNumber}{2}){
						$spotMedianF_B_storage{$peptideNumber}{2} = $spotMedianF_B;
						#print "stored peptide #: $peptideNumber\t2\t$spotMedianF_B_storage{$peptideNumber}{2}\n";
					}
					elsif (!exists $spotMedianF_B_storage{$peptideNumber}{3}){
						$spotMedianF_B_storage{$peptideNumber}{3} = $spotMedianF_B;
						#print "stored peptide #: $peptideNumber\t3\t$spotMedianF_B_storage{$peptideNumber}{3}\n";
					}
					else {
						$spotMedianF_B_storage{$peptideNumber}{4} = $spotMedianF_B;
						#print "stored peptide #: $peptideNumber\t4\t$spotMedianF_B_storage{$peptideNumber}{4}\n";
					}
				}
			}
		}
		else {
			next
		}
	}
	#Normalize the F-B values to the brightest non-control spot on each array
	#print "key1\tkey2\tequation\tvalue\n";
	#print "PeptideID\tSorted1\tSorted2\tSorted3\tSorted4\tMedian\n";
	foreach my $key1 (sort keys %spotMedianF_B_storage){#peptideID
		my $mid;
		my @normalizedSpotMedianValues;
		foreach my $key2 (sort keys %{$spotMedianF_B_storage{$key1}}){#replicateID
			my $normalizedSpotMedian = $spotMedianF_B_storage{$key1}{$key2}/$max_F_B;
			#print "$key1\t$key2\t$spotMedianF_B_storage{$key1}{$key2} / $max_F_B\t$normalizedSpotMedian\n";
			push @normalizedSpotMedianValues,$normalizedSpotMedian;#add each median value to array
			$medianF_B{$sampleName}{$key1} = $medianPeptideValue;
			#print "key1: $key1\tkey2: $key2\tvalue: $spotMedianF_B_storage{$key1}{$key2}\n";
		}
		##calculate median value across all spots for each peptide
			my @sorted_values = sort @normalizedSpotMedianValues;#ranks lowest to highest in array
			$mid = int @sorted_values/2;
			#print "sorted: $sorted_values[0]\t$sorted_values[1]\t$sorted_values[2]\t$sorted_values[3]\tmid: $mid\n";
			if (@sorted_values % 2){
				$medianPeptideValue = $sorted_values[$mid];
			}
			else{
				$medianPeptideValue = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
			}
			$medianF_B{$sampleName}{$key1} = $medianPeptideValue;
			#print "$key1\t$sorted_values[0]\t$sorted_values[1]\t$sorted_values[2]\t$sorted_values[3]\t$medianF_B{$sampleName}{$key1}\n";
	}
}
close (ORIGFILE);

##output all results for serum infected with each virus type to separate file
#generate column headers
my $summaryOutputFile = "Summarized_Peptide_Results.tsv";
open (SUMMARY, ">$summaryOutputFile") || die "$summaryOutputFile: $!\n";
print SUMMARY "Characterized_Result\tSample_ID\t";#to offset for sample labels
#foreach my $MMkey1 (sort keys %medianF_B){#Sample_ID
for (my $i = 0; $i < $sortedLength; $i++){#peptideID
	#print "$i\t";
	print SUMMARY "$sortedPeptides[$i]\t";
	#print "$i\t$sortedPeptides[$i]\n";
}
#}
print SUMMARY "\n";

#generate peptide names in column headers
print SUMMARY "\t\t";
#foreach my $MMMMkey1 (sort keys %medianF_B){#peptide metadata
for (my $j = 0; $j < $sortedLength; $j++){#peptideID
	#print "$j\t";
	print SUMMARY "$name{$sortedPeptides[$j]}\t";
}
#}
print SUMMARY "\n";

#generate peptide sequences in column headers
print SUMMARY "\t\t";
#foreach my $MMMkey1 (sort keys %medianF_B){#peptide metadata
for (my $k = 0; $k < $sortedLength; $k++){#peptideID
	#print "$k\t";
	print SUMMARY "$sequence{$sortedPeptides[$k]}\t";
}
#}
print SUMMARY "\n";

#iterate over sorted hash again to print values to file in same order as column headers
foreach my $Mkey1 (sort keys %medianF_B){#Sample_ID
	print SUMMARY "$characterizedSamples{$Mkey1}\t$Mkey1\t";
	#print "$characterizedSamples{$Mkey1}\t$Mkey1\t";
	for (my $l = 0; $l < $sortedLength; $l++){#peptideID
		if (exists $medianF_B{$Mkey1}{$sortedPeptides[$l]}){
			print SUMMARY "$medianF_B{$Mkey1}{$sortedPeptides[$l]}\t";
			#print "$Mkey1\t$l\t$sortedPeptides[$l]\t$medianF_B{$Mkey1}{$sortedPeptides[$l]}\n";`
		}
		else {
			print SUMMARY "\t";
			my $position1 = $l+1;
			print "Peptide $position1 not present in results file for $Mkey1\n";
		}
	}
	print SUMMARY "\n";
}
close (SUMMARY);

#TSRI Code	Characterized
__END__
H1	DV1
H2	DV1
H3	DV1
H4	DV1
H5	DV1
H6	DV1
H7	DV1
H8	DV2
H9	DV2
H10	DV2
H11	DV2
H12	DV2
H13	DV2
H14	DV2
H15	DV2
H16	DV3
H17	DV3
H18	DV3
H19	DV3
H20	DV3
H21	DV3
H22	Dengue Negative 
H23	Dengue Negative 
H24	ZIKV, DENV
H25	ZIKV, DENV
H26	ZIKV, DENV
H27	ZIKV, DENV
H28	ZIKV, DENV
H29	ZIKV, DENV
H30	ZIKV, DENV
H31	ZIKV, DENV
H32	ZIKV, DENV
L1	DV3
L2	DV2
L3	DV3
L4	DV2
L5	DV2
L6	DV2
L7	DV3
L8	DV3
L9	DV2
L10	DV2
L11	DV2
L12	DV2
L13	DV2
L14	DV2
L15	DV2
L16	DV2
L17	DV3
L18	DV2
L19	DV2
L20	DV2
L21	DV2
L22	DV3
L23	DV2
L24	ZIKV
L25	ZIKV
L26	ZIKV
L27	ZIKV
L28	ZIKV
L29	ZIKV
L30	ZIKV
L31	ZIKV
L32	ZIKV
C1	CHKV
C2	CHKV
C3	CHKV
C4	CHKV, DENV
C5	CHKV, DENV
C6	CHKV, DENV
C7	CHKV, DENV
C8	CHKV, DENV
C9	CHKV, DENV
C10	CHKV, DENV
C11	CHKV, DENV
C12	CHKV 
C13	CHKV, DENV
C14	Unknown
C15	CHKV 
C16	CHKV, DENV
C17	CHKV, DENV
C18	CHKV, DENV
C19	CHKV, DENV
C20	CHKV, DENV
C21	CHKV, DENV
C22	CHKV, DENV
C23	CHKV, DENV
C24	DENV
C25	DENV
C26	CHKV, DENV
C27	CHKV, DENV
C28	CHKV, DENV
C29	CHKV, DENV
C30	CHKV, DENV
C31	CHKV, DENV
C32	CHKV, DENV
C33	CHKV, DENV
C34	CHKV, DENV
C35	CHKV, DENV
C36	CHKV, DENV
C37	CHKV, DENV
C38	CHKV, DENV
C39	CHKV, DENV
C40	CHKV, DENV
C41	DENV
C42	DENV
C43	DENV
C44	DENV
C45	DENV
C46	DENV
C47	DENV
C48	DENV
C49	DENV
C50	DENV
C51	WNV
C52	WNV
C53	WNV
C54	WNV
C55	WNV
C56	WNV
C57	WNV
C58	WNV
C59	WNV
C60	WNV
C61	WNV
C62	WNV
C63	WNV
C64	ZIKV
C65	ZIKV
C66	ZIKV
C67	ZIKV
C68	ZIKV
C69	ZIKV
C70	ZIKV
C71	ZIKV
C72	ZIKV
C73	ZIKV
C74	Unknown
C75	ZIKV
C76	ZIKV
