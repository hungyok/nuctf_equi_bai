# Nucleosome & TF occupancy using SEMs
We present a computational methods to compute nucleosome and TF occupancy using statistical equilibrium models (SEMs). 
## Input data
### Genomic sequence
Download the 16 chromosomes (12071326 bp) of Saccharomyces cerevisiae (S288C) in fasta format from SGD https://www.yeastgenome.org/.
### TF motif database
Use their recommended cutoff for every TF and download in total 104 budding yeast TF position-weight-matrices from http://stormo.wustl.edu/ScerTF/. The list of TFs “listbai_all.txt” has all the 104 TF names (column1) and their sizes (column2). All the 104 PWMs are ordered in tandem (same ordering as in listbai_all.txt), with the four columns ordered as A/C/G/T and saved in a file called “wmsbai_data_all.txt”. The two files listbai_all.txt and wmsbai_data_all.txt can be found in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all).
### Nucleosome Mnase-seq data
transf_yy_data ---> It has code that converts log-scale to linear-scale occupancy 
Lee et al.
## Processing binding energy and NDR 
nuc_energy ---> Nucleosome energy for all 16 chromosomes

tf_energy_all ---> It has all the 16-chromosome sequences (source: SGD), 104 TF PWMs, TF index list, and codes to compute TF energy. The TF energy is dumped in the folder: tf_energy_all\Etf_allmat_chr. Due to large data size, we have reported only the first three chromosome (Etf_chr1.mat, Etf_chr2.mat, and Etf_chr3.mat) and each has only the first 10 TFs from the listbai_all.txt.

ndr_call ---> Codes to compute NDR positions 
## SEM and optimization
NucTF ---> Has codes for optimization and nucleosome/TF occupancy with TFs. It also has examples to compute the final occupancy are present.

NucRemod ---> Has codes for optimization and nucleosome/TF occupancy with TFs and remodelers. It also has examples to compute the final occupancy.
## Output data
Nucleosome occupancy, TF occupancy.










