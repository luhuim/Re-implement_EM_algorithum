#!/usr/bin/python3
"""
Title: alculating motif matrix using EM
Date: 2023 March 11th
Author: Huimin Lu

Description:
    This script is to caculating DNA motif matrix using EM (Expection Maximization)
    An input file is a fasta file that includes couple of unaligned DNA sequences
        and the length of each sequence must be same
    And input parameter is:
        X=0.5
        W: the length of motif

Function:
    


Procedure:
    1, setting the start points for EM, which is initial p-matrix
    2, Iterated E-step, calculating matrix-Z using p-matrix
    3, Iterated M-step, calculating matrix-p using matrix-Z
    4, finally got matrix-Z that converged
    5, write input file for tool "weblogo"
Usage:
    python main.py dataset/Dataset3.fasta output_file/dataset_3.txt -kmer 6 
"""

#Function that will be used
def initial_matrix_p (sub_seq, X=0.5, A=0.25,T=0.25,C=0.25,G=0.25):
    
    ##pandas is required
    ## sub_seq is the subsequence that is randomly chosen
    ## X is the probability of motif, default if 0.5
    ## this raw index of matrix is "A T C G"
    alphabet={"A":'0',
              "T":'1',
              "C":'2',
              "G":'3'}
    data={}
    #adding first column into dictrionary "data"
    data[str(0)]=[A,T,C,G]
    for i in range(0,len(sub_seq)):
        col=[0,0,0,0]
        letter=sub_seq[i]#one letter in subsequence
        col[int(alphabet[letter])]=X#The index where X should loacate
        col=[round((1-X)/3,2) if x!=X else X for x in col  ]
        # print(col)
        data[str(i+1)]=col
    # print("data",data)
    # Creates pandas DataFrame.
    P0 = pd.DataFrame(data, index=['A',
                                   'T',
                                   'C',
                                   'G'])
    # print(P0)
    return P0




"""
EM(DNA_sequences,last_P, kmer, A,T, C, G )
    This function is to calculate next matrix_P on the top of last matrix_P
    "DNA_sequences" is a list that includes all seqeunces in dataset
    "last_P" is a dataframe, matrix_p that was generated last time
    "kmer" is the length of motif
    "A","T","C","G" are background probability
"""
def EM(DNA_sequences,last_P, kmer, A,T, C, G ):
    print("inside the function")
    # #Step-E
    alphabet={"A":'0',
              "T":'1',
              "C":'2',
              "G":'3'}
    #create an empty Z-matrix
    l=len(DNA_sequences[0])#"l" is the length of one DNA seqeunce
    Z_col=[str(i) for i in range(1,l-kmer+2)]#column: 1,2,3,4,5....
    # print("Z_col",Z_col)
    # Z_raw=["seq"+str(i) for i in range(1,s+1)]#seq1,seq2,seq3
    # Z={}#key is sequance name called "seq1".... value is all Z values in this row
    matrix_Z = pd.DataFrame(columns=Z_col)#create an empty dataframe
    
    j=[i for i in range(0,l-kmer+1)]#the starting index of motif in sequences
    # print("j",j)
    h=1
    for seq in sequences:#look are one DNA seqeunce one by one
        print(seq)
        zi=[]#"zi" include one row of Z-value
        for i in j:#scan every subsequence in one seqeunce
            #one subsequence, calculating one cumprod() and add into Z-matrix
            zii=[]#multipy all element in this list,result is one z value
            # current_subseq=seq[i:i+W]
            for x in range(0,len(seq)):#scan every letter in one sequence, "x" is index 
                
                if x<i or x>=i+kmer:#if this letter is in front of motif
                    current_letter=seq[x]#current letter
                    pi=last_P["0"][current_letter]#one pi value
                    zii.append(pi)#add pi into list "zi"
                else:#if this letter is motif
                    current_letter=seq[x]#current letter
                    pi=last_P[str(x-i+1)][int(alphabet[current_letter])]#one pi value
                    #df.iloc
                    # print(type(last_P))
                    # print(last_P)
                    # print("current index",x)
                    # print("index in kmaer",x-i+1)
                    # print("current_letter",current_letter)
                    # print("corresponding column",last_P[str(x-i+1)])
                    # print("corresponding element",pi)
                    zii.append(pi)#add pi into list "zi"
            # print(zii)
            result = np.prod(np.array(zii))#cumprod the value, "result" is one Z-value in Z_matrix
            # print("multiple,one element",result)
            # print("current_subseq",zii)
            zi.append(result)#"zi" is one row of zi value in Z_matrix
            zii=[]
        # Z["seq"+str(h)]=zi#add one row in dictionary "Z"
        # print("one line of Z matrix",zi)
        # print(Z_col)
        #add one line into Z matrix
        row_name="seq"+str(h)
        normalized_zi = preprocessing.normalize([zi])#normalization z value in one line
        # print("normalized row in Z matrix",normalized_zi)
        # print(normalized_zi)# "normalized_zi" has double square bracket, it looks like [[ value]]
        matrix_Z.loc[row_name] = normalized_zi[0]#add normalized z value into matrix Z
        # matrix_Z = pd.DataFrame(np.array([zi]), columns=Z_col,index=[row_name])
        # print(matrix_Z)
        h+=1
    print("***************matrix-Z this round********************")
    print(matrix_Z)
##M-step
#create a matrix_P
    matrix_P = pd.DataFrame(columns = [str(i) for i in range (1,kmer+1)],
            index = ['A', 'T', 'C','G'])

    for c in range(0,kmer):#"c" is index inside of motif(column in P-matrix)
        for base in ["A","T","C","G"]:#scan "base" one by one(row in P-matrix)
            #print(base)
            pii=[]
            for a in range (0,len(DNA_sequences)):#each sequence
                line=DNA_sequences[a]#"line" is one sequence in dataset
                for b in range(0,l-kmer+1):# 
                    subsequence=line[b:b+kmer]#"subsequence" is motif
                    if subsequence[c]==base:#whether this letter in motif is equal to "base"
                        #add "Z value" into list "pii"
                        pii.append(matrix_Z[str(b+1)]["seq"+str(a+1)])  # row: seq, col: index of that letter b+c 
            # print(pii)
            pi=(sum(pii)+1)/(matrix_Z.values.sum()+4)
            matrix_P[str(c+1)][str(base)]=pi
            pii=[]
    #adding one column for background probability
    matrix_P.insert(0, "0", [A,T,C,G],True)
    print("----------------------matrix P this round---------------------------")
    print(matrix_P)
    return(matrix_P)
#%%
import random
import pandas as pd
import numpy as np
from sklearn import preprocessing
import argparse
import math

#call a parser
parser = argparse.ArgumentParser(description='Replication EM alorgrithm in DNA motif discovery')
#adding parameter
#parser.add_argument("parg") 
parser.add_argument('input_file', type=str,
                    help='input fasta file, it must be standard fasta file')
parser.add_argument('output_file', type=str,
                    help='output txt file, this file will be used to make motif logo')
parser.add_argument('-kmer', type=int,
                    help='the length of motif')
parser.add_argument("-A", type=float,default=0.25,
                    help="background(non motif probability)")
parser.add_argument("-T", type=float,default=0.25,
                    help="background(non motif probability)")
parser.add_argument("-C", type=float,default=0.25,
                    help="background(non motif probability)")
parser.add_argument("-G", type=float,default=0.25,
                    help="background(non motif)probability")
parser.add_argument("-X", type=float,default=0.5,
                    help="probalility of motif among dataset, used in calculate initial p matrix")

# loading parameters
args = parser.parse_args()
 
#call parameter
#print(args.parg)
input_fa=args.input_file
output_txt=args.output_file
#print("echo ={0}".format(args.digit))
W=args.kmer
A=args.A
T=args.T
C=args.C
G=args.G
X=args.X
# print(args.parg)

##import input file, and making subsequences
input_seq=open(input_fa,'r')#input file, W=6

sequences=[]#store original DNA seqeunces
subseqeunces=[]
for line in input_seq:
    if line.startswith(">"):
        pass
    elif len(line.strip())==0:#avoid empty line
        pass
    else:#this line is a DNA seqeunce
        #remove new line sign, transform letter into upper-case and add into list "sequence"
        sequences.append(line.strip().upper())
        l=len(line)#l is length of input DNA seqeunce
j=[i for i in range(0,l-W+1)]#the starting index of motif in sequences
for seq in sequences:
    for i in j:
        subseqeunces.append(seq[i:i+W])
subseqeunces=list(set(subseqeunces))
chosen_subsequence=random.choice(subseqeunces)#randomly choose one subsequence from list "subsequences""
# print(subseqeunces)
#using function "initial_matrix_p" to calculate inintial_P matrix
initial_P=initial_matrix_p (chosen_subsequence)#"initial_P" is first matrix-P
# print(initial_P)
#using function "EM" to calculate first P_matrix
first_matrix_P=EM(sequences,initial_P, W, A,T,C, G )
# print("first_matrix_P",first_matrix_P)
pre=initial_P #P-matrix after first iterate
minus=[]
for i in range(0,40):
    post=EM(sequences,pre, W, A,T,C, G  )
    # print("P-matrix",post)
    diff=post-pre#subtruct two matrix
    # print("diff",diff)
    #make every element squared
    diff_square=diff*diff
    # print("diff_square",diff_square)
    diff_square=np.matrix(diff_square)
    diff_square_sum=diff_square.sum()
    # print("diff_square_sum",diff_square_sum)
    difference=math.sqrt(diff_square_sum)
    # print("difference",difference)
    if difference <0.0001:
        break
    else:
        pre=post
        i+=1



#Writing final P-matrix in output file
post_t = post.T 
output=open(output_txt,"w")
output.write("NA motif0\nXX\nID motif0\nXX\n")
output.write("P0\tA\tT\tC\tG\n")
for i in range(0,W+1):
    line=post_t.iloc[i]
    line=list(line)
    line=[str(round(i,9)) for i in line]
    line='\t'.join(line)
    output.write("P"+str(i+1)+"\t"+line+"\n")
output.write("XX\n//")
output.close()
input_seq.close()
# #%%
# import numpy as np
# matrix = np.array([[1, 2, 3, 4],
#     [5, 6, 7, 8],
#     [9, 10, 11, 12]])
# print(matrix*matrix)



