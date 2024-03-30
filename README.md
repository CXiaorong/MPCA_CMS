# MPCA_CMS
In this project, a Multi-core Parallel Co-evolution Algorithm is proposed for solving the CMS model (MPCA_CMS). The input is a binary mutation matrix A, a connected PPI network Q and a parameter K. The output is a set of genes corresponding to submatrix M. In the MPCA_CMS method, due to the large number of genes and samples, it takes much time in calculating individual fitness. Therefore, the calculation of individual fitness is partitioned into a set of independent computation tasks and assigned to multi-core processors to execute simultaneously.

## Operating environment：
Windows 10，R3.4.1

## Input datas: a weighted non-binary mutation matrix A, a PPI network Q, a parameter K; <br>
* binary mutation matrix: data\GBM_SNVdata_440.csv ;    data\OVCA_SNVdata_2547.csv ;  data\THCA_SNVdata_3420.csv <br>
Example of A input to algorithm,  Their rows represent the same set of cancer samples, and their columns represent two sets of genes.<br>
![image](https://github.com/CXiaorong/MPCA_CMS/assets/105973069/4aadb1e2-15eb-4b46-bc66-2757d3a21348)

* network matrix: data\GBM_network_440.csv ;   data\OVCA_network_2547.csv ;  data\THCA_SNVdata_3420.csv <br>
Example of Q file input to algorithm, Both their rows and columns represent genes.
![image](https://user-images.githubusercontent.com/105973069/169654040-765489f4-7d48-44f3-89d3-73e204380797.png)

## Output: a set of genes corresponding to submatrix M;	
	A gene set and its corresponding fitness function value.
 For example: GBM dataset, k=6, the result is: "CDKN2A" "EGFR"   "ERBB2"  "MDM2"   "PIK3R1" "TP53" , "2.31534391534392"
 
 ## Steps:
### 1.Install and load the parallel package.
    install.packages('doParallel')
    library(doParallel)
### 2.First select the functions and run it.
![image](https://github.com/CXiaorong/MPCA_CMS/assets/105973069/eb96a0b3-e91d-452e-a6d3-f301be06c424)

### 3.load data
![image](https://github.com/CXiaorong/MPCA_CMS/assets/105973069/eee4559d-ab4d-4917-bd2f-46a7a85792ff)

### 4.Initialization parameters.
![image](https://github.com/CXiaorong/MPCA_CMS/assets/105973069/7ad30ecc-33e8-4652-963c-e66eff5df032)

### 5.Iterative loop
![image](https://github.com/CXiaorong/MPCA_CMS/assets/105973069/d00bf869-cd98-47bf-84f2-0b9e5e3ef63b)


 
