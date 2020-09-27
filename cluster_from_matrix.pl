#!/usr/bin/env perl

##############################################################
##  Author: Jason Caravas                                   ##
##  Description: Uses DBSCAN algorithm to naively cluster   ##
##                distance or similarity matrices.          ##
##############################################################

use strict;
use warnings;

use Math::Complex;
use File::Path qw/make_path/;
use Getopt::Long;
use Cwd;

my $buckets_default = 1;

GetOptions (
	'matrix=s'			=>	\my $matrix_file,						#Required input matrix
	'prefix=s'			=>	\my $prefix,							#Required file prefix
	'epsilon:f'			=>	\my $epsilon,							#Optional epsilon
	'buckets:i'			=>	\(my $buckets = $buckets_default),		#Optional buckets
	'minpts:i'			=>	\my $minpts, 							#Optional MinPts
	'distance'			=>	\my $distance,							#Matrix is already distance
	'scale:i'			=>	\(my $scale = 100), 					#Optional Scale	 - 100 for percent, 1 for fraction
	'rowcol:s'			=>	\(my $rowcol = "row"),					#Where are matrix labels?  "row", "col", or "both"
	'outdir:s'			=>	\my $outdir,							#Where to store results files
	'rlib:s'			=>	\my $rlib,								#Where is local R library for fpc and dbscan?
	'help'				=>	\my $help,								#Run help and exit		
) or die "Invalid options passed to $0\n";

my $cwd = getcwd();

if (defined $help) {
	Usage();
}
CheckParams ();
PrintParams ();

#  Parse matrix and convert to distance matrix if necessary
my $dist_matrix = ReadMatrix ($matrix_file, $rowcol, $distance, $scale);

if (defined $outdir) {
	unless (-d $outdir) {
		make_path $outdir;
	}
	chdir $outdir or die $!;
}

#  Print out matrix
my $dist_matrix_file = $prefix . "_dist_matrix.txt";
unlink $dist_matrix_file;
PrintMatrix ($dist_matrix_file, $dist_matrix);
my $returned_files;
push (@$returned_files, ["Reformatted distance matrix", $dist_matrix_file]);

#  Calculate $minpts
$minpts = CalcMinPts ($minpts);

#  Calculate $epsilon
$epsilon = CalcEpsilon ($dist_matrix_file, $prefix, $epsilon, $minpts, $buckets);

#  Generate clusters and figures
my $calc_cluster_files = CalcClustersDbscan ($dist_matrix_file, $prefix, $epsilon, $minpts, $buckets); 
push (@$returned_files, @$calc_cluster_files);

#  Make a readable cluster file - pass dist matrix and cluster file produced by dbscan
my @translate_clusters_returns = TranslateClusters ($dist_matrix, $returned_files -> [-1] -> [1]);

push (@$returned_files, $translate_clusters_returns[0]);
my $cluster_count =  $translate_clusters_returns[1];

print "\n$cluster_count clusters + noise (\"cluster 0\") identified using minimum cluster size of $minpts and clustering distance of $epsilon\n"; 

print "\nFiles produced:\n";
foreach my $file (@$returned_files) {
	print "\t", $file -> [0], ":  ", $file -> [1], "\n";
}
print "\n";

sub CheckParams {
	unless (defined $matrix_file) {
		print "Input matrix must be specified with --matrix argument\n";
		Usage();
	}
	unless (defined $prefix) {
		print "File prefix must be specified with --prefix argument\n";
		Usage();
	}
	unless (-f $matrix_file) {
		print "Specified matrix file $matrix_file is not a valid file name\n";
		Usage();
	}
	if ((defined $minpts) && (uc $minpts eq "AUTO")) {
		$minpts = '';
	}	
	if ((defined $epsilon) && (uc $epsilon eq "AUTO")) {
		$epsilon = '';
	}
	if ((defined $epsilon) && ($epsilon <= 0)) {
		print "Specified epsilon value must be greater than 0\n";
		Usage();
	}
	if ((defined $buckets) && ($buckets < 1)) {
		print "Specified buckets value must be greater than or equal to 1\n";
		Usage();
	}
	if ((defined $minpts) && ($minpts <= 0)) {
		print "Specified minpts value must be greater than 0\n";
		Usage();
	}
	unless (($scale == 1) || ($scale == 100)) {
		print "Specified scale value must be exactly 1 or 100\n";
		Usage();
	}
	$rowcol = lc $rowcol;
	unless (($rowcol eq "row") || ($rowcol eq "col") || ($rowcol eq "both")) {
		print "Specified rowcol value must be \"row\", \"col\", or \"both\"\n";
		Usage();
	}			
}

sub PrintParams {
	print "Provided arguments:\n";
	print "\tmatrix: $matrix_file\n";
	print "\tprefix: $prefix\n";
	if ($buckets != $buckets_default) { 
		print "\tbuckets: $buckets (Overridden)\n";
	}
	else { 
		print "\tbuckets: $buckets\n";
	}
	if (defined $epsilon) { 
		print "\tepsilon: $epsilon (Overridden)\n";
	}
	else { 
		print "\tepsilon: Calculate\n";
	}
	if (defined $minpts) { 
		print "\tminpts: $minpts (Overridden)\n";
	}
	else { 
		print "\tminpts: Calculate\n";
	}
	if (defined $distance) { 
		print "\tdistance: Input matrix is a distance matrix (Overridden)\n";
	}
	else { 
		print "\tdistance: Input matrix is a similarity matrix\n";
	}
	if ($scale == 100) {
		print "\tscale: Value scaling $scale\n";
	}
	elsif ($scale == 1) {
		print "\tscale: Value scaling $scale (Overridden)\n";
	}
	if ($rowcol eq "row") {
		print "\trowcol: Matrix labels are in \"$rowcol\"\n";
	}
	else {
		print "\trowcol: Matrix labels are in \"$rowcol\" (Overridden)\n";
	}
	if (defined $outdir) {
		print "\toutdir: Output directory set to $outdir\n";
	}
	else {
		print "\toutdir: Output directory defaulting to $cwd\n";
	}	
	if (defined $rlib) {
		print "\trlib: Local R library location set to $rlib\n";
	}
	else {
		print "\trlib: Local R library location not set\n";
	}	
	
	
	print "\n";		
							
}	


sub ReadMatrix {
	my $file = shift;
	my $rowcol = shift;
	my $distance = shift;
	my $scale = shift;
	
	my $matrix;
	my $labels;
	
	print "Parsing matrix input file . . . ";

	if ($rowcol eq "col") {
		my $matrix_info = ParseColMatrix ($file);
		$matrix = $matrix_info -> {"matrix"};
		$labels = $matrix_info -> {"labels"};
	}
	elsif ($rowcol eq "row") { 
		my $matrix_info = ParseRowMatrix ($file, "0");
		$matrix = $matrix_info -> {"matrix"};
		$labels = $matrix_info -> {"labels"};
	}
	elsif ($rowcol eq "both") { 
		my $matrix_info = ParseRowMatrix ($file, "1");
		$matrix = $matrix_info -> {"matrix"};
		$labels = $matrix_info -> {"labels"};
	}	
			
	my $return_matrix;
	
# 	my $prec_100 = Math::BigFloat -> new (100);
	
	##  Walk through matrix and fill in square matrix
	for (my $i = 0; $i < @$matrix; $i++) {
		
		##  Fill in diagonal (0 distance)
		$return_matrix -> [$i] -> [$i] = 0;		
		
		for (my $j = 0; $j < @{$matrix -> [$i]}; $j++) {
			my $val = $matrix -> [$i] -> [$j];
			#  Unless specified as distance matrix, assume similarity matrix		
			unless (defined $distance) {
				$val = $scale - $val;
			}			
			$val = SmallNums ($val);
 			
			##  Fill in matrix
			$return_matrix -> [$i] -> [$j] = $val;
			$return_matrix -> [$j] -> [$i] = $val;
		}
	}
	##  Unshift labels onto top row
	unshift (@$return_matrix, $labels);
	print "Done!\n";
	
	return $return_matrix;
}

sub ParseRowMatrix {
	my $file = shift;
	my $skip_header = shift;
	
	my $matrix;
	my $labels;
	
	open (my $in_fh, "<", $file) or die;
	
	#  If row and col headers in file, just throw away col headers
	if ($skip_header) {
		my $header = <$in_fh>;
	}
			
	while (my $line = <$in_fh>) {
		$line =~ s/\R+//g;
		chomp $line;
		my @cols = split (/\s+/, $line);
		my $id = shift @cols;
		push (@$labels, $id);
		map {$_ = SmallNums($_)} @cols;	
		push (@$matrix, \@cols);
	}
	my $return_val;
	$return_val -> {"matrix"} = $matrix;
	$return_val -> {"labels"} = $labels;
	close $in_fh;

	return $return_val;	
}
sub ParseColMatrix {
	my $file = shift;
	
	my $matrix;
	my $labels;

	open (my $in_fh, "<", $file) or die;
	my $header = <$in_fh>;
	$header =~ s/\R+//g;
	chomp $header;
	$header =~ s/^\s+//g;
	my @header_labels = split (/\s+/, $header);
	$labels = \@header_labels;	
			
	while (my $line = <$in_fh>) {
		$line =~ s/\R+//g;
		chomp $line;
		my @cols = split (/\s+/, $line);
		map {$_ = SmallNums($_)} @cols;
		push (@$matrix, \@cols);
	}
	my $return_val;
	$return_val -> {"matrix"} = $matrix;
	$return_val -> {"labels"} = $labels;
	close $in_fh;
		
	return $return_val;	
}

sub PrintMatrix {
	my $file = shift;
	my $matrix = shift;
	
	print "Printing distance matrix to $file . . . ";	
	open (my $out_fh, ">", $file) or die $!;
	foreach my $row (@$matrix) {
		print $out_fh join ("\t", @$row), "\n";
	}
	print "Done!\n";
	close $out_fh;
}

sub CalcMinPts {
	my $minpts = shift;
	
	print "Calculating minpts . . . ";
	unless (defined $minpts) {
		#  Get count of elements
		my $count = @$dist_matrix;
		#  Subtract column header
		$count -= 1;
		my $unrounded = log ($count);
		$minpts = int ($unrounded);
		
		#  Round up
		if ( $unrounded - $minpts > 0) {
			$minpts++;
		}
		print "Value $minpts (rounded up from $unrounded)\n";
		return $minpts;
	}
	else {
		print "Overridden at command line - $minpts\n";
		return $minpts;
	}
}

sub CalcEpsilon {
	my $dist_matrix_file = shift;
	my $prefix = shift;
	my $epsilon = shift;
	my $minpts = shift;
	my $buckets = shift;
	
	print "Calculating epsilon . . . ";
	unless (defined $epsilon) {
		my $rscript_file = $prefix . "_calc_epsilon.rscript";
		open (my $out_fh, ">", $rscript_file) or die $!;
		
		if (defined $rlib) {
			print $out_fh "library (\"dbscan\", lib.loc\"$rlib\")\n";
		}
		else {
			print $out_fh "library (\"dbscan\")\n";
		}
		print $out_fh "matrix <- read.table(file=\"$dist_matrix_file\",header=TRUE,sep=\"\\t\")\n";
		print $out_fh "dist_matrix <- as.dist(matrix)\n";
		print $out_fh "dist_matrix <- as.dist(matrix)\n";		
		print $out_fh "knndist <- dbscan::kNNdist(dist_matrix, $minpts, bucketSize=$buckets, splitRule=\"SUGGEST\")\n";
# 		print $out_fh "knndist\n";
		print $out_fh "write.table(knndist, sep=\"\\t\", row.names=TRUE)\n";
		
		close $out_fh;
		my $return = `Rscript $rscript_file`;
		my @rows = split (/\R+/, $return);
		my @knn_list;
		foreach my $row (@rows) {
			#  Skip "x" row
			unless ($row =~ m/\d+/) {
				next;
			}
			$row =~ s/\R+//g;
			chomp $row;
			my @cols = split (/\s+/, $row);
			push (@knn_list, $cols[-1]);
		}
		my $epsilon = FindApex(\@knn_list);
		print "Value $epsilon\n";
		return $epsilon;
				
	}
	else {
		print "Overridden at command line - $epsilon\n";
		return $epsilon;
	}
}

sub FindApex {
	my $kNN_list = shift;
	my @kNN_list = @$kNN_list;
	@kNN_list = reverse sort {$a <=> $b} @kNN_list;
# 	print join (",", @kNN_list), "\n";
	
	my $length = @kNN_list;
	my $b = $kNN_list[0];
	my $m = ($kNN_list[-1] - $kNN_list[0])/$length;
	my $tm = -$m;
	

	
# 	print "\n\tSlope = ", sprintf("%.12f", $m), "\n\tY-intercept = $b\n";
# 	print "\tTangent slope = ", sprintf("%.12f", $tm), "\n"; 
	
	my $best_y;
	my $best_dist = 0;
	for (my $i = 0; $i < @kNN_list; $i++) {
		 my $val = $kNN_list[$i];
		# Calculate X intercept of tangent
		my $tb = -(($tm * $i) - $val);
# 		print "Tangent intercept for $val = $tb\n";
		
		# Solve for x position of intercept
		# A little algebra and we have x = ($tb - $b)/(2 * $m)
		my $x_int = ($tb - $b)/(2 * $m);

		# Solve for y position of intercept
		# Just fill in y=mx+b equation
		my $y_int = ($m * $x_int) + $b;

#		Sanity Check :-)				
# 		my $y2 = ($tm * $x_int) + $tb;
# 		print "$y_int = $y2\n";		
		
		# Calc distance between curve and line
		# Use known curve location ($i, $val) as 1st value
		# Use intercept ($x_int, $y_int) as 2nd value
		
# 		no bignum;
		my $dist = sqrt ( (($x_int - $i) ** 2) - (($y_int - $val) ** 2) );
# 		print "$dist\n";
		if ($dist > $best_dist) {
			$best_dist = $dist;
			$best_y = $val;
# 			print "\t$best_y\n";
		}	
		
	}
	return $best_y;
	
}
sub CalcClustersDbscan {
	my $dist_matrix_file = shift;
	my $prefix = shift;
	my $epsilon = shift;
	my $minpts = shift;
	my $buckets = shift;
	
	my @returned_files;
	
	print "Calculating clusters . . . ";
	my $rscript_file = $prefix . "_calc_clusters.rscript";
	open (my $out_fh, ">", $rscript_file) or die $!;
	
	if (defined $rlib) {
		print $out_fh "library (\"dbscan\", lib.loc\"$rlib\")\n";
	}
	else {
		print $out_fh "library (\"dbscan\")\n";
	}
	
	print $out_fh "matrix <- read.table(file=\"$dist_matrix_file\",header=TRUE,sep=\"\\t\")\n";
	print $out_fh "dist_matrix <- as.dist(matrix)\n";
	
	#  knn table
	my $knn_table_file = $prefix . "_knn_table.txt";
	unlink $knn_table_file;
	my $returned_file = ["k-Nearest Neighbor table", $knn_table_file];
	push (@returned_files, $returned_file);  
	
	print $out_fh "knndist <- dbscan::kNNdist(dist_matrix, $minpts, bucketSize=$buckets, splitRule=\"SUGGEST\")\n";
	print $out_fh "write.table(knndist, file=\"$knn_table_file\", sep=\"\\t\", row.names=TRUE)\n";
	
	#knn distance plot
	my $knn_dist_plot = $prefix . "_knn_dist_plot.pdf";
	unlink $knn_dist_plot;	
	$returned_file = ["k-Nearest Neighbor distance plot", $knn_dist_plot];
	push (@returned_files, $returned_file);  
	
	print $out_fh "pdf(file=\"$knn_dist_plot\")\n";
	print $out_fh "dbscan::kNNdistplot(dist_matrix, $minpts, bucketSize=$buckets, splitRule=\"SUGGEST\")\n";
	print $out_fh "abline(h=$epsilon, col=\"red\")\n";
 	print $out_fh "title(main=\"Epsilon=$epsilon\")\n";	
	print $out_fh "dev.off()\n";
	
	#Switch from dbscan to fpc
	print $out_fh "detach(\"package:dbscan\", unload=TRUE)\n";
	if (defined $rlib) {
		print $out_fh "library (\"fpc\", lib.loc\"$rlib\")\n";
	}
	else {
		print $out_fh "library (\"fpc\")\n";
	}
	
	#dbscan plot
	my $dbscan_plot = $prefix . "_dbscan_plot.pdf";
	unlink $dbscan_plot;	
	$returned_file = ["DBSCAN plot", $dbscan_plot];
	push (@returned_files, $returned_file); 	

	print $out_fh "pdf(file=\"$dbscan_plot\")\n";
	print $out_fh "dbscan_clus <- fpc::dbscan(dist_matrix, $epsilon, MinPts=$minpts, method=\"dist\", showplot=1)\n";
	print $out_fh "dev.off()\n";
	
	#dbscan cluster file
	my $dbscan_clusters = $prefix . "_dbscan_clusters.txt";
	unlink $dbscan_clusters;	
	$returned_file = ["DBSCAN clusters", $dbscan_clusters];
	push (@returned_files, $returned_file); 	

	print $out_fh "write.table(dbscan_clus\$cluster, file=\"$dbscan_clusters\",sep=\"\\t\", row.names=TRUE)\n";
			
	`Rscript $rscript_file`;
	
	print "Done!\n";
	return \@returned_files;
}			

sub TranslateClusters {
	my $dist_matrix = shift;
	my $cluster_file = shift;
	
	my $labels = $dist_matrix -> [0];

	my $final_clusters = $prefix . "_final_clusters.txt";
	unlink $final_clusters;	
	my @returned_file = ("Final clusters", $final_clusters);
	
	#  Yes, kludging this in here is bad form, but I didn't want to do a 2nd pass to find it
	my $max_cluster = 0;	
	
	open (my $in_fh, "<", $cluster_file) or die $!;
	open (my $out_fh, ">", $final_clusters) or die $!;
	
	while (my $line = <$in_fh>) {
		$line =~ s/\R+//g;
		chomp $line;
		
		if ($line =~ m/\"(\d+)\"\t(\d+)/) {
			my $index = $1 - 1;	#  Convert 1 based R table to 0 based array index
			my $cluster = $2;
			print $out_fh $labels -> [$index], "\t$cluster\n";
			if ($cluster > $max_cluster) {
				$max_cluster = $cluster;
			}
		}
	}
	
	close $in_fh;
	close $out_fh;
	
	return (\@returned_file, $max_cluster);
}

sub SmallNums {
	#  Forcing small nums into 10 digit decimal notation	
	my $val = shift;
	if ($val < 1) {
		$val = sprintf("%.10f", $val);
	}
	return $val;
}

sub Usage {
	print "\nUSAGE:\n";
	print "\nperl cluster_from_matrix.pl --matrix <MATRIX FILE> --prefix <PREFIX> [options]\n";
	print "\n";
	print "This script takes a similarity or distance matrix and clusters the results\n";
	print "using the \"Density-based spatial clustering of applications with noise\"\n";
	print "(DBSCAN) algorithm.\n";
	print "\n";
	print "All clustering parameters are automatically determined using commonly accepted\n";
	print " estimators,unless overridden by the user.\n";
	print "\n";
	print "This software requires the \"dbscan\" and \"fpc\" R packages to be installed.\n";
	print "\n";
	
	print "PARAMETERS:\n";
	printf '  %-10s  %s', "--matrix", "User provided similarity or distance matrix [REQUIRED]";
	print "\n\n";
	printf '  %-10s  %s', "--prefix", "Output file naming prefix [REQUIRED]";
	print "\n\n";
	printf '  %-10s  %s', "--help", "Display this help text";
	print "\n\n";	
	printf '  %-10s  %s', "--outdir", "Location to save output files. [DEFAULT: Current directory]";
	print "\n\n";
	printf '  %-10s  %s', "--rlib", "Location of local R libraries. [DEFAULT: No path necessary]";
	print "\n\n";
	printf '  %-10s  %s', "--distance", "Provided matrix is already a distance matrix. <TOGGLE> [DEFAULT: OFF]";
	print "\n\n";
	printf '  %-10s  %s', "--rowcol", "Distance matrix headings are located on \"row\", \"col\", or both.";
	print "\n";
	printf '  %-10s  %s', "", "[DEFAULT: \"row\"]";
	print "\n\n";							
	printf '  %-10s  %s', "--scale", "If a similarity matrix is provided, the scale is used to calculate the distance";
	print "\n";
	printf '  %-10s  %s', "",        "matrix (distance = scale - similarity value).If the similarity matrix is expressed";
	print "\n";
	printf '  %-10s  %s', "",        "in terms of percentages, \"100\" should be specified. If expressed in decimal";
	print "\n";
	printf '  %-10s  %s', "",        "proportions, \"1\" should be used. This parameter is ignored if \"--distance\" is set.";
	print "\n";
	printf '  %-10s  %s', "",        "[DEFAULT: 100]";
	print "\n\n";	
	printf '  %-10s  %s', "--minpts", "Minimum dbscan cluster size expressed as an integer. [DEFAULT: AUTO]";
	print "\n\n";
	printf '  %-10s  %s', "--epsilon", "Epsilon density value for dbscan clustering [DEFAULT: AUTO]";
	print "\n\n";
	printf '  %-10s  %s', "--buckets", "Bucket size for dbscan clustering [DEFAULT: 1]";
	print "\n\n";	
		
	exit;
}
	
	