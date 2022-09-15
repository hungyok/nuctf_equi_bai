# Nucleosome & TF occupancy using SEMs
We present a computational method to compute nucleosome and transcription factor (TF) occupancies using statistical equilibrium models (SEMs). For complete details on the use and execution of this method, please refer to [Kharerin & Bai, 2021](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008560). 
## Input data
### Genomic sequence
Download the 16 chromosomes (12071326 bp) of *Saccharomyces cerevisiae* (S288C, genome version R64-3-1) in fasta format from SGD https://www.yeastgenome.org/. The downloaded genomic sequence can be found [here](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all/sgd_genome).
### TF motif database
File "[listbai_all.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all/listbai_all.txt)" contains the name and motif size of 104 TFs that are found in the budding yeast nuclei in normal YPD growth condition ([Yan, Chen & Bai](https://doi.org/10.1016/j.molcel.2018.06.017)). File "[wmsbai_data_all.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all/wmsbai_data_all.txt)" shows the combined PWMs of these 104 TFs in tandem (same ordering as in listbai_all.txt), with the four columns representing the propensity of finding A/C/G/T at each motif position. The recommended score cutoff is shown in "[PWM_TF_cutoff.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all)". These files can be found in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all).
### Nucleosome occupancy data
The original nucleosome occupancy data “[analyzed_data_complete_bw20.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/transf_yy_data)” was provided in Lee et al.’s supplementary data ([Lee, et al., 2007](https://www.nature.com/articles/ng2117)). The first four sequences are non-genomic DNA and excluded them from the dataset. We rename the new file as "[raw_lee_chr1_16.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/transf_yy_data)" which has chr1 to chr16. This data, which are in log scale, were further normalized, linearized, and saved as "[yy1_lee.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call/dataFolder)". The normalization and linearization can be conducted with the following code:
```
> expY = importdata('yourpath/transf_yy_data/raw_lee_chr1_16.txt');
> [c,γ] = transf_yy_data(expY);
> d_lee=importdata('yourpath/transf_yy_data/skimpped_data.txt');
> d_indx=importdata('yourpath/transf_yy_data/skimpped_data_index.txt');
> [x1_lee,y1_lee] = linear_occup(expY,c,γ,d_lee,d_indx);
> save yy1_lee.mat x1_lee y1_lee;
```
The MATLAB code transf_yy_data.m and linear_occup.m can be found in the folder [transf_yy_data](https://github.com/hungyok/nuctf_equi_bai/tree/main/transf_yy_data). The final occupancy "yy1_lee.mat" can be found in the folder [dataFolder](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call/dataFolder).

## Processing binding energy and NDR 
### Nucleosome energy
The sequence-dependent binding energy for nucleosome is obtained using the software developed by Kaplan & Segal ([Kaplan, et al., 2009](https://www.nature.com/articles/nature07667)), which can be downloaded from the following link:

 https://github.com/KaplanLab/NucleosomeModel
 
 The link has all the instructions on how to use the software. Briefly, download the code “NucleosomeModel-master.zip” on Linux 64-bit machine. Enter the following commands on the terminal to generate the energy files: 
 ```
> make install
> perl nucleosome_prediction.pl -raw_binding -t example -s input.fa -p raw_output -tab
```
Where input.fa is the DNA sequence in fasta format. One sequence per file, e.g., chr1.fa has sequence for chromosome 1 only. To be consistent with our codes, we suggest all the sequences be in caps. The raw_output.tab is the output energy (log-score) file in tab format (can easily be open by any text editor like "notepad"). This energy is stored as [E_Em.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/nuc_energy) in the folder [nuc_energy](https://github.com/hungyok/nuctf_equi_bai/tree/main/nuc_energy). 
### TF energy
The binding energy of a TF is obtained by scanning its PWM along the genomic sequences (both forward and reverse-complement strands) and converting into position-dependent energy. The maximum of the two possible energies at a position is the TF energy.

The code “[tf_binding_pot.m](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all)” uses the genome sequence, wmsbai_data_all.txt, and listbai_all.txt to obtain the energy profiles. Run the commands below to get all the 104 TF energies:
```
> path1 = 'yourpath/tf_energy_all/sgd_genome/'; % enter full path
> path2 = 'yourpath/tf_energy_all/Etf_allmat_chr/'; % enter full path
> [Emtf,k] = tf_binding_pot(path1,path2);
```
The genomic sequence (e.g., sgd_chr1.fa for chromosome 1 and so on) from the folder “sgd_genome” is fed using the input parameter “path1”. The code then spews out TF energy profiles for all the chromosomes (Etf_chr1.mat to Etf_chr16.mat) and the average TF energies (Emtfall.mat), which are stored in the folder “Etf_allmat_chr” indicated by “path2”. Note that we have reported here only three energy files, Etf_chr1.mat to Etf_chr3.mat, due to large file sizes! Each energy file (e.g., Etf_chr1.mat for chromosome 1) is a matrix file where it contains the binding energy for each TF at different chromosomal locations. The code, the subfolders, and the input/output files can be found in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all). 
### NDR annotation
To call and locate nucleosome-depleted-regions (NDRs) from the genome-wide nucleosome occupancy data, we carry out the following two steps:

**ndr_cut.m**: Using “yy1_lee.mat” as input data, this program identifies regions below certain occupancy threshold. In our current approach, we used eight different thresholds: threshold(i)=0.8-i×SD where i=1,2,…,8, and SD is the standard-deviation calculated from genome-wide nucleosome occupancy (yy1_lee.mat), which yield eight output files: NDR_1.mat, NDR_2.mat, …, NDR_8.mat. Set SD=0.0678 (default) and run the following code with input i=1,2,…,8.
```
> load yourpath/ndr_call/dataFolder/yy1_lee.mat;
> % step above generates two variables, x1_lee and y1_lee
> NDR_i = ndr_cut(x1_lee, y1_lee, i); 
```
**ndr_pos_cal2.m**: This program compares the boundaries in NDR_1.mat, NDR_2.mat, etc., picks the NDRs with steep edges, merges nearby NDRs, and force the nucleosome occupancies in NDRs to be zero. Copy all the NDR_i.mat into the "dataFolder":
```
> path1 ='yourpath/ndr_call/dataFolder/';
> [yy3A, ndrpos_chrA] = ndr_pos_cal2(path1);
```
The final output includes the NDR positions ([ndrpos_chrA.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call)) and the modified nucleosome occupancy ([yy3A.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call)). The occupancy, yy3A.mat, is renormalized to 80% and is recorded as [yy3A_lee.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call). Use this modified occupancy for future model optimization. These codes and the input/output data can be found in the folder [ndr_call](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call).
## SEM and optimization
The SEM partition function is full of this type of term: p<sub>i</sub>=c<sub>t/N</sub>e<sup>-γ<sub>t/N</sub>E<sup>i</sup><sub>t/N</sub></sup> where p<sub>i</sub> is an individual particle (nucleosome, N, or TF, t) Boltzmann weight. We use a modified Nelder-Mead simplex algorithm with Simulated Annealing to optimize the scaling factors c<sub>t/N</sub> and γ<sub>t/N</sub> by calculating the the root-square-mean deviation (RMSD) of nucleosome occupancy between the expiremental data (yy3A_lee.mat) and our model.

### NucTF
To evaluate individual TF contribution to NDRs, we first optimize (c,γ) for an individual TF when only the concerned TF is present in the model. We use a modified Nelder-Mead simplex algorithm with Simulated Annealing to optimize these parameters.

#### Open the folder “[simplexM_tf1](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_tf1)” to find:
1. Input and output files that are saved in folder “input” and “output” respectively.
2. A log file “log.txt” to record the optimization process.
3. simplex_SA_singleTF.m, the main simplex engine that executes the parameter optimization.
4. occupxfunc.m, a sub-function that computes nucleosome occupancy and RMSD.
5. sortf.m, a sub-function that sorts RMSDs
#### Edit directory paths:
```
> path1 ='yourpath/nuc_energy/'; 
> path2 = 'yourpath/tf_energy_all/';    
> path3 ='yourpath/NucTF/simplexM_tf1/'
 ```
#### Run the optimization code:
```
> simplex_SA_singleTF(tfx,path1,path2,path3);
```
#### Input parameters:
1. The "tfx" is the interest of TF index listed in "lisbai_all.txt". It is any number from 1 to 104.  
2. The "path1" loads nucleosome energy “E_Em.mat” located in the folder [nuc_energy](https://github.com/hungyok/nuctf_equi_bai/tree/main/nuc_energy).
3. The "path2" loads TF energy files Etf_chr1.mat to Etf_chr16.mat, and Emtfall.mat located in the folder [Etf_allmat_chr](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all/Etf_allmat_chr). It also loads listbai_all.txt, a list of 104 TF names and motif sizes, located in the folder [tf_energy_all](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all). 
4. The "path3" loads yy3A_lee.mat and rand_genome.mat, where the former is the reference nucleosome occupancy data, and the latter consists of two elements. The first element is the number of chromosomes that accounts for 70% of the genome (we tune the parameters on 70% of the genome and use the rest 30% as an independent test). The second element is a random 1D array of the 16 chromosomes. Both yy3A_lee.mat and rand_genome.mat are in the folder simplexM_tf1/[input](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_tf1/input). The path3 also loads tf_xhIIa.mat (simplexM_tf1/input/[initialp](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_tf1/input/initialp)), a vector of eight numbers: an initial guessed values of (c<sub>N</sub>;γ<sub>N</sub>;c<sub>t</sub>;γ<sub>t</sub>) (the first four numbers) and corresponding small deviations about the latter values (the last four numbers). 
#### Output files:
1. simplex_xval.txt stores sorted (c,γ).
2. simplex_fval.txt stores sorted RMSDs.
3. simplex_tval.txt stores latest values of difference between the best and worst RMSDs, temperature, and RMSD.
4. simplex_Temp.txt reports the content of “simplex_tval.txt” per simplex step.
	
The file simplex_xval.txt is a 4x5 matrix of c and γ. The four rows of the matrix are c<sub>N</sub>, γ<sub>N</sub>, c<sub>t</sub>, and γ<sub>t</sub>, respectively. The first column of the matrix is our optimized parameters. These files are dumped in the folder [output](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_tf1/output).

After carrying out the above one-TF optimization for all the TFs, we rank the 104 TFs based their contribution to the NDR probability (P<sub>NDR</sub>) by running the "[occupx.m](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile)" and store the ranking in another file called "[tfindx.txt](https://github.com/hungyok/nuctf_equi_bai/tree/main/tf_energy_all)". Below, we describe a model where we incorporated the top 30 TFs. Accordingly, we have in total of 62 unknown parameters (31 pairs of (c,γ)'s), of which 60 is for TFs and two for nucleosome.

#### Open the folder “[simplexM_top30](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_top30)” to find:
1. Input and output parameters that are saved in folder "input" and "output" respectively. The folder "output" contains similar files as in simplexM_tf1. We also put some optimized parameters into "output_optimized" for comparison purpose. The "output" is the running folder where every output is dumped here. The output for real or trial run, abrupt kill, or mistakes by users are dumped here. 
2. A log file “log.txt” to record the optimization process.
3. simplex_SA_multiTF.m, the main simplex engine that executes the parameter optimization.
4. occupxfunc.m, a sub-function that computes nucleosome occupancy and RMSD.
5. sortf.m, a sub-function that sorts RMSDs.

#### Edit directory paths:
```
> path1 ='yourpath/nuc_energy/'; 
> path2 = 'yourpath/tf_energy_all/';    
> path3 ='yourpath/NucTF/simplexM_top30/'; 
```
#### Run the optimization code:
```
> simplex_SA_multiTF(tfx,path1,path2,path3) ;
```
#### Input parameters:
Here we have four input parameters: tfx, path1, path2, and path3. The functions of these inputs are the same as described in the previous section, except that tfx here is a 1D array of TF indices. For example, to consider top 30 TFs that contribute to NDRs, the tfx should be the first 30 numbers listed in the "tfindx.txt" file. Path3 loads "top30_xhIIa.mat" which is a 1D array containing the initial guessed values of (c,γ)s of nucleosome and the 30 TFs. Here, the first 62 elements are the (c,γ)s and the second 62 are the small displacements needed to build simplex vertices. This file can be found in the folder [initialp](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_top30/input/initialp). 
#### Output files: 
1. simplex_xval.txt stores sorted (c,γ).
2. simplex_fval.txt stores sorted RMSDs.
3. simplex_tval.txt stores latest values of difference between the best and worst RMSDs, temperature, and RMSD.
4. simplex_Temp.txt reports the content of “simplex_tval.txt” per simplex step.

Similar file description as in simplexM_tf1 except that the file simplex_xval.txt is a 62x63 matrix of c and γ with the first and second rows of the matrix are c<sub>N</sub> and γ<sub>N</sub> respectively, and the rows 3–32 are c<sub>t</sub>’s and rows 33–62 are γ<sub>t</sub>’s. The first column of the matrix is our optimized parameters and will be used to determine the final genome-wide occupancy profile (see below). 

Number of search steps at a given temperature “M” in simplex_SA.m can be changed to adjust the optimization quality/time. Currently, at each of the five “temperature” (which controls the parameter search range), we gradually decrease the number of steps to be M=100, 83, 67 ,50, and 33. Larger M may lead to better convergence but cost more computational time. For the current prescribed set of M (=100, 83, 67 ,50, 33), we repeated step 3 nine times and each time starting from the last output (Figure 3). One simplex step takes &sim; 125 seconds. We recommend users to run a cumulative of at least &sim; 3000 simplex steps &sim; 4 days of simulations. 

Open the folder [occup_profile](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile) and run the code to compute the occupancy profiles:
```
> path1 ='yourpath/nuc_energy/'; 
> path2 = 'yourpath/tf_energy_all/';    
> path3 ='yourpath/ndr_call/';
> load yourpath/NucTF/simplexM_top30/output_optimized/simplex_xval.txt;
> xfit = simplex_xval (:,1); TF = 1; Eseq = 1;
> load yourpath/NucTF/tf_energy_all/tfindx.txt;
> tfx = tfindx(1:30);
> [O] = occupx(TF,Eseq,xfit,tfx,path1,path2,path3);
> save O_good.mat O; % save O_bad.mat O; when loading simplex_xval.txt from folder “output”
```
The logical parameters “TF” and “Eseq” can take on the value of 0 or 1. For example, TF=1 means that TF information is considered in the model otherwise TF=0, and Eseq=1 means that nucleosome sequence-specific energy is considered and otherwise set Eseq=0. The output is an occupancy data, O.mat, comprising of 16x1 cell arrays (one for each chromosome), which includes the occupancies for nucleosome and each TF considered in the model. The code “occupx.m” is preloaded with nucleosome energy, TF energy, and TF index lists (listbai_all.txt and tfindx.txt). The output and figure of two example cases when "TF=0, Eseq=0" and "TF=0, Eseq=1" can be found in the folders “[example1](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example1)” and “[example2](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example2)”, respectively. The output of a model example with TF=1 (top 30 TFs) and Eseq=1 can be found in the folder “[example3](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example3)”. All the codes and input/output files and folders can be found in the folder [simplexM_top30](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/simplexM_top30).

### NucRemod
This version of the model considers both the effect from TF binding and nucleosome remodeling (NR). We assume that the action of NR is to deform the nucleosome energy landscape near the TF-bound sites in a TF occupancy-dependent manner. NR is present only when TF occupancy > tfcut (set tfcut=0.0022 computed as an genome-wide average of occupancy for all TFs in the model NucTF). We assume that an NR modifies the energy landscape as a Gaussian deformation with height “h” in kT/0.212 and width “w” in bp (free parameters). In the presence of multiple TFs adjacent to each other, we add up individual deformation energies. 
#### Open the folder “[simplexM_remod](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucRemod/simplexM_remod)” to find: 
1. An “input” folder containing a file “pos_octf.mat” that records positions and occupancy of TFs when the TF occupancy cutoff is tfcut=0.0022. It is generated by the code “occup_tfs.m” by running the following steps:
```
> path1 ='yourpath\nuc_energy\'; 
> path2 = 'yourpath\tf_energy_all\';   
> tfcutx = 0.0022;
> load yourpath/NucRemod/simplexM_remod/input/avgxval_top30.txt; 
> xfit=avgxval_top30.txt; 
> load yourpath/NucTF/tf_energy_all/tfindx.txt;
> tfx = tfindx(1:30);
> [pos_octf] = occup_tfs(xfit,tfx,tfcutx,path1,path2);
> save pos_octf.mat pos_octf;
```
The format of pos_octf.mat is a 16x3 cell matrix. Each row in the cell represents a different chromosome (first row is chr1, 2nd row is chr2, and so on). First column represents TF binding position, while the 2nd and 3rd column represents the occupancy and identity of the TF. This file is used by occupR.m, occupRfunc.m, and tf_cluster.m. 

2. Two output folders for fitting parameters: "output" and "output_optimized". Parameters generated during the fitting process will be stored in "output". While "output_optimized" contains some precomputed best fitting parameters. 
3. An output log file “log.txt” that records the optimization process.
4. simplex_SA_remod.m, the simplex engine with Simulated Annealing.
5. occupRfunc.m, a sub-function that calls tf_cluster.m and occup_nucs.m to construct remodeling deformation energy, calculate nucleosome occupancy, and compute RMSD. This step is repeated nx = 10 times to allow combinatorial binding of multiple TFs (see step below). The final nucleosome occupancy represents an average between these configurations. Larger nx can be used to achieve better averaging. 
6. tf_cluster.m, a sub-function that gathers the binding configurations of TFs using the information provided by "pos_octf.mat". For isolated TFs, the binding status is determined by comparing their occupancy (from model NucTF) to a random number. When multiple TFs bind adjacent to each other, we allow maximally two TFs to overlap, and their individual binding status are determined by the same random number generator. The output of this function was used to generate an altered energy landscape in occupRfunc.m.
7. occup_nucs.m, a sub-function that computes nucleosome occupancy (one chromosome at a time) using the altered energy landscape. 
8. sortf.m, a sub-function that sorts a set of RMSDs.
#### Edit directory paths:
```
> path1 ='yourpath\nuc_energy\'; 
> path2 = 'yourpath\tf_energy_all\';    
> path3 ='yourpath\NucTF\simplexM_top30\';
> path4 ='yourpath\NucRemod\simplexM_remod\';
```
#### Run the optimization code:
```
> simplex_SA_remod(tfx,path1,path2,path3,path4);
```
Note that tfx is the same index array as in "simplexM_top30" in model NucTF. The prescribed set of M used here is M=12,12,10,10,6. We repeated step 3 four times and each time starting from the last output. One simplex step takes &sim; 1 hour. We recommend users to run a cumulative of at least &sim; 200 simplex steps &sim; 8 days of simulations. The input and output files are same as in NucTF, except here the file simplex_xval.txt is a 4x5 matrix. The first column of the matrix (c<sub>N</sub>, γ<sub>N</sub>, h, and w) is our optimized parameters and the TF parameters (c<sub>t</sub>, γ<sub>t</sub>) obtained in NucTF are used to calculate the occupancy profiles. The occupancy calculation is not as straight forward as in NucTF, but it depends on “nx” to account for the TF-dependent remodeling. Here we set nx=100. Open the folder [occup_profile](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucRemod/occup_profile) and run the code to compute the occupancy profiles:
```
> path1 ='yourpath\nuc_energy\'; 
> path2 = 'yourpath\tf_energy_all\';    
> path3 ='yourpath\NucRemod\simplexM_remod\';
> path4 ='yourpath\ndr_call\';
> load yourpath/NucRemod/simplexM_remod/output_optimized/simplex_xval.txt;
> remod = simplex_xval(:,1); % load the best-fit remodeling parameters
> load yourpath/NucTF/simplexM_top30/output_optimized/simplex_xval.txt;
> xfit = simplex_xval(:,1); % load the best-fit TF parameters
> load yourpath/NucTF/tf_energy_all/tfindx.txt;
> tfx = tfindx(1:30);
> [O] = occupR(remod,xfit,tfx,path1,path2,path3,path4); 
```
The output file “O.mat” is a 16x1 cell. Here we report only nucleosome occupancy per cell and can be found in the folder “[occup_profile](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucRemod/occup_profile)”. The folder also contains a figure showing comparison between a model with TF only and a model with TF+remodeling, and a plotting script. All the codes and input/output files and folders can be found in the folder [simplexM_remod](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucRemod/simplexM_remod). 
## Output data
The main output files from this exercise are the occupancy profiles of nucleosome and TF, and a list of annotated NDRs "[ndrpos_chrA.mat](https://github.com/hungyok/nuctf_equi_bai/tree/main/ndr_call)". The output files "O.mat" without remodeling effect (NucTF) can be found in the folders [example1](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example1), [example2](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example2), and [example3](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucTF/occup_profile/example3). The output file "O.mat" with remodeling effect (NucRemod) can be found in [occup_profile](https://github.com/hungyok/nuctf_equi_bai/tree/main/NucRemod/occup_profile).

 
