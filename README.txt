Created: 23/11/2020
Author: Joanna J. Ilska

Script calculating H matrix from genotypes provided in Plink.raw format, and ordered pedigree. 
The G portion of the matrix calculated according to VanRaden's 1st method. 
Missing genotypes replaced with mean genotype. 

The script takes two arguments:
1) path to pedigree file
2) path to genotypes in Plink .raw format

The IDs are recoded to numeric. 

Produces 3 output files:
<filename>_H.giv - the inverse of H matrix
<filename>_pedRecoded.csv
recode_IDs.csv

Folder includes an example: ped.txt and gen.txt

 The calculation of A matrix is the most time consuming part. For a pedigree with 6,200 individuals, this took about 6h. 
 Could be improved if Ainverse could be read in from an external source.