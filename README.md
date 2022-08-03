# Nucleosome & TF occupancy using SEMs
We present a computational methods to compute nucleosome and TF occupancy using statistical equilibrium models (SEMs). 
## Input data
### Genomic sequence
Download the 16 chromosomes (12071326 bp) of Saccharomyces cerevisiae (S288C) in fasta format from SGD https://www.yeastgenome.org/.
### TF motif database
Use their recommended cutoff for every TF and download in total 104 budding yeast TF position-weight-matrices from http://stormo.wustl.edu/ScerTF/. The list of TFs “listbai_all.txt” has all the 104 TF names (column1) and their sizes (column2). All the 104 PWMs are ordered in tandem (same ordering as in listbai_all.txt), with the four columns ordered as A/C/G/T and saved in a file called “wmsbai_data_all.txt”. The two files listbai_all.txt and wmsbai_data_all.txt can be found in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all).
### Nucleosome Mnase-seq data
Use the processed dataset “analyzed_data_complete_bw20.txt” provided in Lee et al.’s supplementary data ([Lee, et al., 2007](https://www.nature.com/articles/ng2117)). We have also provided a copy in the folder [transf_yy_data](https://github.com/hungyok/nuctf_equi_bai/tree/main/transf_yy_data). The first four sequences are non-genomic DNA and so we excluded them from the dataset. We rename the new file as "raw_lee_chr1_16.txt".

All data points >1.26 are considered outliers which is only 0.0011% of the total data (3017253) and are discarded. Then, the raw nucleosome occupancy data above (S) in log-scale is converted to linear-scale occupancy (Y) using the formula: Y=ce<sup>γS</sup>. We varied c and γ values to satisfy the following conditions: 1) the value of Y at each genomic index should range from 0 to 1, 2) the average 〈Y〉 is 80% (experimental data show that yeast genome is ~80% occupied by nucleosomes in yeast), and 3) γ is maximized to enhance the contrast between high and low nucleosome occupancy. 

The MATLAB code transf_yy_data.m can execute the above step and can be found in the folder [transf_yy_data](https://github.com/hungyok/nuctf_equi_bai/tree/main/transf_yy_data). 
```
> S2 = importdata('raw_lee_chr1_16.txt');
> [c,g] = transf_yy_data(S2);
 ```
The code takes the raw data (raw_lee_chr1_16.txt) as an input and finds the best (c, γ) that satisfy the above conditions. In this example c = 0.8244 and γ= 0.1581.

## Processing binding energy and NDR 
### Nucleosome energy
Sequence-dependent binding energy for nucleosome is obtained using the software developed by Kaplan & Segal ([Kaplan, et al., 2009](https://www.nature.com/articles/nature07667)), which can be downloaded from the following link:

 https://github.com/KaplanLab/NucleosomeModel
 
 The link has all the instructions on how to use the software. Briefly, download the code “NucleosomeModel-master.zip” on Linux 64-bit machine. Enter the following commands to generate energy files: 
 ```
> make install
> perl nucleosome_prediction.pl -raw_binding -t example -s input.fa -p raw_output -tab
```
Where input.fa is the DNA sequence in fasta format. One sequence per file. E.g., chr1.fa has sequence for chromosome 1 only. To be consistent with our codes, we suggest all the sequences be in caps. The raw_output.tab is the output energy (log-score) file in tab format (can easily be open by any text editor like "notepad"). 

nuc_energy ---> Nucleosome energy for all 16 chromosomes

tf_energy_all ---> It has all the 16-chromosome sequences (source: SGD), 104 TF PWMs, TF index list, and codes to compute TF energy. The TF energy is dumped in the folder: tf_energy_all\Etf_allmat_chr. Due to large data size, we have reported only the first three chromosome (Etf_chr1.mat, Etf_chr2.mat, and Etf_chr3.mat) and each has only the first 10 TFs from the listbai_all.txt.

ndr_call ---> Codes to compute NDR positions 
## SEM and optimization
NucTF ---> Has codes for optimization and nucleosome/TF occupancy with TFs. It also has examples to compute the final occupancy are present.

NucRemod ---> Has codes for optimization and nucleosome/TF occupancy with TFs and remodelers. It also has examples to compute the final occupancy.
## Output data
Nucleosome occupancy, TF occupancy.










