## Nucleosome occupancy using SEM

## The second largest heading


The repository folder "nuctf_equi_bai" has the data and codes to compute nucleosome abd TF occupancy. It has the following subfolders:

transf_yy_data ---> It has code that converts log-scale to linear-scale occupancy 

ndr_call ---> Codes to compute NDR positions  

nuc_energy ---> Nucleosome energy for all 16 chromosomes

tf_energy_all ---> It has all the 16-chromosome sequences (source: SGD), 104 TF PWMs, TF index list, and codes to compute TF energy. The TF energy is dumped in the folder: tf_energy_all\Etf_allmat_chr. Due to large data size, we have reported only the first three chromosome (Etf_chr1.mat, Etf_chr2.mat, and Etf_chr3.mat) and each has only the first 10 TFs from the listbai_all.txt.

NucTF ---> Has codes for optimization and nucleosome/TF occupancy with TFs. It also has examples to compute the final occupancy are present.

NucRemod ---> Has codes for optimization and nucleosome/TF occupancy with TFs and remodelers. It also has examples to compute the final occupancy.
