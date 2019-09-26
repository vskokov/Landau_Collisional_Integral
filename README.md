# Landau Collisional Integral

## Dependencies: 

GSL 

## Compilation 
g++ -o fp_exp.x -O3 FP_expansion.cpp -lgsl -lgslcblas 

## Run 
fp_exp.x 

## Plot Results
xmgrace -nxy output_data.csv

The first column is t; the second is T; the remaining are c_1, ..., c_Lmax

## Sample output 

L_max = 10: output_data_10.csv

L_max = 20: output_data_20.csv
