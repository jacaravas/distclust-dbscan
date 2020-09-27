I have provided an example alignment and distance matrix for the NCBI ortholog group PCLA_564588 (adenine deaminase), https://www.ncbi.nlm.nih.gov/proteinclusters/564588  While this example is for distance calculated from an amino acid alignment, the method works on any sort of distance or similarity matrix formatted several ways.

Test command:
perl cluster_from_matrix.pl –matrix test/PCLA_564588_distance.txt –prefix PCLA_564588 –distance –rowcol col
