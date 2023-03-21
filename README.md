# Find_motif_in_dataset
Write software that finds motifs in a dataset (generate your own data). See an example of an existing application https://meme-suite.org/meme/tools/meme 
## Version
```
0.0.1
```
# Create environment and install package
Create a conda environment called `binp29`
```
conda create --name binp29 python=3.9
```
Package used in this project
```
conda   22.9.0
weblogo 3.7
python  3.9
pandas  1.4
numpy   1.21
```
Install `weblogo` in conda environment
```
conda install -c conda-forge weblogo
```
# Usage
Create a directory called `find_motif`, enter into this directory and clone this repository inside this directory.
```
mkdir find_motif
cd find_motif
git clone https://github.com/luhuim/find_motif.git
```
## Runing the program
create position weight matrix
```
python src/main.py data/Dataset3.fasta results/dataset_3.txt -kmer 6
python src/main.py data/Dataset6.fasta results/dataset_6.txt -kmer 6 
python src/main.py data/Dataset7.fasta results/dataset_7.txt -kmer 6 
```
create motif logo using `weblogo`
```
weblogo -D table -F 'png' -A 'dna' -a 'ATCG' -s large -c classic <results/dataset_3.txt> results/logo_3.png || exit
weblogo -D table -F 'png' -A 'dna' -a 'ATCG' -s large -c classic <results/dataset_6.txt> results/logo_6.png || exit
weblogo -D table -F 'png' -A 'dna' -a 'ATCG' -s large -c classic <results/dataset_7.txt> results/logo_7.png || exit
```


