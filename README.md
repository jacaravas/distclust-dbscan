# distclust-dbscan
Naive clustering of distance and similarity matrices. Supports several formats of square or diagonal matrix with row or column labels. Optional auto-calculation of epsilon using a "knee"-finding approach.

Dependencies:
Perl (Tested with 5.22.1)
	Math::Complex
	File::Path
	Getopt::Long
	Cwd
R (Tested with 3.6)
	dbscan library
	fpc library


In-program help:
USAGE:

perl cluster_from_matrix.pl --matrix <MATRIX FILE> --prefix <PREFIX> [options]

This script takes a similarity or distance matrix and clusters the results
using the "Density-based spatial clustering of applications with noise"
(DBSCAN) algorithm.

All clustering parameters are automatically determined using commonly accepted
 estimators,unless overridden by the user.

This software requires the "dbscan" and "fpc" R packages to be installed.

PARAMETERS:
  --matrix    User provided similarity or distance matrix [REQUIRED]

  --prefix    Output file naming prefix [REQUIRED]

  --help      Display this help text

  --outdir    Location to save output files. [DEFAULT: Current directory]

  --rlib      Location of local R libraries. [DEFAULT: No path necessary]

  --distance  Provided matrix is already a distance matrix. <TOGGLE> [DEFAULT: OFF]

  --rowcol    Distance matrix headings are located on "row", "col", or both.
              [DEFAULT: "row"]

  --scale     If a similarity matrix is provided, the scale is used to calculate the distance
              matrix (distance = scale - similarity value).If the similarity matrix is expressed
              in terms of percentages, "100" should be specified. If expressed in decimal
              proportions, "1" should be used. This parameter is ignored if "--distance" is set.
              [DEFAULT: 100]

  --minpts    Minimum dbscan cluster size expressed as an integer. [DEFAULT: AUTO]

  --epsilon   Epsilon density value for dbscan clustering [DEFAULT: AUTO]

  --buckets   Bucket size for dbscan clustering [DEFAULT: 1]

