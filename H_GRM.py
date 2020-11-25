## Author: Joanna J. Ilska
## Date: 19/11/2020

"""
Script calculating an H GRM matrix from genotypes in a Plink .raw format, and a pedigree (csv format with a header).
For the G portion, VanRaden 1 is used.  

Pedigree needs to be presorted!

Produces inverse of H in <filename>_H.giv
Recoded pedigree in filename_pedRecoded.csv
and recode_IDs.csv

The IDs in grm are recoded, starting from 1. The script produces a file <filename>_grm_recodedIDs.txt where the first column
is the new IDs included in the GRM and the second column is the original sample IDs. 
"""

import numpy as np
import pandas as pd

import sys as sys
import os as os

def calc_NRM(ped):
	## Create an empty dataframe of the appropriate size
	
	zeros= np.zeros((len(ped),len(ped)))
	np.fill_diagonal(zeros,1)
	A = pd.DataFrame(data = zeros, columns=ped['id'],index=ped['id'])
	
	for i in ped['id']:
		## Extract sire and dam ID for individual i
		sire_i = ped.loc[ped['id']==i]['sire'].values[0]
		dam_i = ped.loc[ped['id']==i]['dam'].values[0]        
		for j in [x for x in ped['id'].values if x>=i]:
			## Extract sire and dam ID for individual j
			sire_j = ped.loc[ped['id']==j]['sire'].values[0]
			dam_j = ped.loc[ped['id']==j]['dam'].values[0]
			if sire_j!=0 and dam_j==0:
				if i!=j:
					stmp = A.loc[i,sire_j]
					A.loc[i,j]=0.5*stmp
					A.loc[j,i]=A.loc[i,j]
			elif sire_j==0 and dam_j!=0:
				if i!=j:
					dtmp = A.loc[i,dam_j]
					A.loc[i,j]=0.5*dtmp
					A.loc[j,i]=A.loc[i,j]
			elif sire_j!=0 and dam_j!=0:
				if i==j:
					jtmp = A.loc[sire_j,dam_j]
					A.loc[i,j]=1+0.5*jtmp
				else:
					ttmp = A.loc[i,sire_j]+A.loc[i, dam_j]
					A.loc[i,j]=0.5*ttmp
					A.loc[j,i]=A.loc[i,j]
			else:
				pass
	return A
	
def VanRadenI(loci):
	'''
	Function calculating GRM using VanRaden's first method, where the product of ZZ' is divided by sum of 2pq, 
	i.e. scaling the G to be analogous to A. 
	Takes as input a numpy array where rows are samples and columns are SNPs.
	'''
	
	## Remove non-informative SNPs
	m1 = np.isnan(loci)
	m2 = loci[0]==loci       
	loci = loci[:,~(m1|m2).all(0)]
	
	## Calculate mean genotype
	mean = np.nanmean(loci, axis=0)
	mean[np.isnan(mean)] = 1
	mean[mean == 0] = 0.00001
	
	## Read through the genotypes are replace missing values with the mean genotype
	for r in range(0, loci.shape[1]):
		loci[np.isnan(loci[:, r]),r] = mean[r]
	
	## Converting mean genotype to matrix shape
	mv = np.repeat([mean],loci.shape[0] , axis=0)
	
	## Calculate the centred genotypes
	W = np.subtract(loci, mv)
	
	## Calculate the sum of 2pq
	freq = mean/2
	pq2=sum(2*freq*(1-freq))
	
	## Calculate the G matrix
	G = np.matmul(W, W.T)
	G = G/pq2
	
	return G
    
if __name__ == "__main__":
	if len(sys.argv) == 3:
		ped_file = sys.argv[1]
		gen = sys.argv[2]
		
		print("Calculating H matrix for individuals in pedigree: ", ped_file)
		print("using genotypes from: ", gen)
		
		## Read in the pedigree
		ped = pd.read_csv(ped_file)
		ped.columns = ['id','sire','dam']
		
		## Remove any potential entries where id=0
		ped = ped.loc[ped['id']!=0]
		ped = ped.loc[ped['id']!='0']
		
		
		## Recode the IDs
		rec_ID = {}
		for i, name in enumerate(ped['id'].values):
			rec_ID[name]=i+1
		
		ped['new_id']=ped['id'].map(rec_ID)
		ped['new_sire']=ped['sire'].map(rec_ID)
		ped['new_dam']=ped['dam'].map(rec_ID)
		ped = ped.drop(columns=['id','sire','dam'])
		ped.columns = ['id','sire','dam']
		ped = ped.fillna(0)
		ped = ped.astype('int64')
		
		
		## Calculate the Ainverse for the whole pedigree
		A = calc_NRM(ped)
		A = A.astype('float')
		
		hinv = np.linalg.inv(np.array(A))
		hinv = pd.DataFrame(data=hinv, index=ped['id'], columns=ped['id'])
		
		
		## Read in the genotypes
		df = pd.read_table(gen, delim_whitespace=True, na_values='NA')
		df = df.replace('K2902B', 'K2902')
		
		sub22=[rec_ID[x.replace('r','').replace('rr','')] for x in df['IID'].values]
		
		## Drop the unnecessary columns
		loci = np.array(df.drop(columns=['FID','IID','PAT','MAT','SEX','PHENOTYPE']))
		print(loci.shape)      
		
		## Calculate G  
		G = VanRadenI(loci)
		## Add a constant to a diagonal
		np.fill_diagonal(G, G.diagonal() + 0.001)
		
		## Calculate the scaling factors
		X = np.array([[1,np.mean(G.diagonal())],[1,np.mean(G)]])
		y=np.array([np.mean(np.array(A.loc[sub22,sub22]).diagonal()), np.mean(np.array(A.loc[sub22,sub22]))])
		s = np.matmul(np.linalg.inv(X),y)
		
		## Scale the G
		Gscaled=s[0]+s[1]*G
		
		## Obtain the inverse of g
		ginv = np.linalg.inv(np.array(Gscaled))
		ginv = pd.DataFrame(data=ginv, index=sub22, columns=sub22)
		
		## Obtain inverse of A22
		ainv22 = np.linalg.inv(np.array(A.loc[sub22,sub22]))
		ainv22 = pd.DataFrame(data=ainv22, index=sub22, columns=sub22)
		
		## Calculate the scaled H22
		df2 = hinv.loc[sub22,sub22]+ginv.loc[sub22,sub22]-ainv22.loc[sub22,sub22]
		
		## Update the matrix
		hinv.update(df2)
		
		## Turn into lower triangular, keep only non zero entries
		stacked = hinv.mask(np.triu(np.ones(hinv.shape),1).astype(bool)).stack().to_frame()
		stacked.index.names=['id1','id2']
		stacked.reset_index(level='id1')
		stacked = stacked.reset_index()
		stacked.rename(columns={0:'Hinv'}, inplace=True)
		
		stacks = stacked.loc[stacked['Hinv']!=0]
		
		outname=os.path.basename(gen).split('.')[0]
		stacks.to_csv("{}_H.giv".format(outname), sep="\t", index=False, header=False)
		
		ped.to_csv("{}_pedRecoded.csv".format(outname), index=False)
		
		w = open("recode_IDs.csv", "w")
		for key, val in rec_ID.items():
			w.write('{},{}\n'.format(key, val))
		
		w.close()
				
	else:
		print("Provide the name of the input .raw file and which GRM to calculate (1 or 2)")
	
    