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
Where input.fa is the DNA sequence in fasta format. One sequence per file, e.g., chr1.fa has sequence for chromosome 1 only. To be consistent with our codes, we suggest all the sequences be in caps. The raw_output.tab is the output energy (log-score) file in tab format (can easily be open by any text editor like "notepad"). 
### TF energy
The binding energy of a TF is obtained by scanning the PWM along the genomic sequences (both forward and reverse-complement strands) and converting into position-dependent energy. The maximum of the two possible energies at a position is the TF energy.

The code “tf_binding_pot.m” uses the genome sequence, wmsbai_data_all.txt, and listbai_all.txt to obtain the energy profiles. Run the commands below to get all the 104 TF energies:
```
> path1 = '/tf_energy_all/sgd_genome/'; % enter full path
> path2 = '/tf_energy_all/Etf_allmat_chr/'; % enter full path
> [Emtf,k] = tf_binding_pot(path1,path2);
```
The genomic sequence (e.g., sgd_chr1.fa for chromosome 1 and so on) from the folder “sgd_genome” is fed using the input parameter “path1”. The code then spews out TF energy profiles for all the chromosomes (Etf_chr1.mat to Etf_chr16.mat) and the average TF energies (Emtfall.mat), which are stored in the folder “Etf_allmat_chr” indicated by “path2”. Each energy file (e.g., Etf_chr1.mat for chromosome 1) is a matrix file where it contains the binding energy for each TF at different chromosomal locations (see the sample data below for four TFs on a few locations on chr I).  The code, the subfolders, and the input/output files can be found in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all). 
```
(Genomic             (TF)
index)	       ABF1	 CBF1	  MCM1	  RAP1
100000       -14.95	-24.07	-22.59	-17.99
100001	      -22.04	-13.19	-21.58	-17.71
100002       -21.29	-13.93	-22.54	-22.67
100003	      -15.56	-23.67	-19.63	-19.17
```
### NDR annotation
To call and locate nucleosome-depleted-region (NDR) we require two steps:

ndr_cut.m: Using “yy1_lee.mat” as input data, this program identifies regions below certain occupancy threshold. In our current approach, we used eight different thresholds: threshold(i)=0.8-i×SD where i=1,2,…,8, and SDis the standard-deviation calculated from genome-wide nucleosome occupancy (yy1_lee.mat), which yield eight output files: NDR_1.mat, NDR_2.mat, …, NDR_8.mat. Set SD=0.0678 (default) and run the following code with input i=1,2,…,8. Note that yy1_lee.mat is directly loaded in the code “ndr_cut.m”. If the code fails to run, specify the right directory path to “load” commands.

ndr_pos_cal2.m: Using “yy1_lee.mat, NDR_1.mat, NDR_2.mat, …, NDR_8.mat” as input (directly loaded into the program), this program generates the final NDR positions (ndrpos_chrA.mat) and the modified nucleosome occupancy map (yy3A.mat). ndrpos_chrA.mat records the start and end index of each NDR on all 16 chromosomes, and yy3A.mat has the same format as yy1_lee.mat.

## SEM and optimization
NucTF ---> Has codes for optimization and nucleosome/TF occupancy with TFs. It also has examples to compute the final occupancy are present.

NucRemod ---> Has codes for optimization and nucleosome/TF occupancy with TFs and remodelers. It also has examples to compute the final occupancy.
## Output data
Nucleosome occupancy, TF occupancy.










