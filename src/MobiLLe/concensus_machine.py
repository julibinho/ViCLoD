"""
Copyright (c) 2019 Bishnu Sarker (bishnukuet@gmail.com), Nika Abdollahi, Juliana Silva Bernardes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""



# This file contains functions to compute the consensus sequence 
# using vector based computations.
# For queries please contact: bishnukuet@gmail.com
import numpy as np


AA2NUM={"A":0,  "D":2,   "C":1  , "E":3 ,  "F":4,   "G":5,   "H":6,   "I":7,   "K":8,"L":9,
        "M":10,   "N":11,   "P":12 ,  "Q":13,   "R":14,   "S":15,   "T":16,   "V":17,   "W":18,   "Y":19, "*" :20, ".":21,"X":22}
NUM2AA={0:"A",1:"C",2:"D",3:"E",4:"F",5:"G",6:"H",7:"I",8:"K",9:"L",
        10:"M",11:"N",12:"P",13:"Q",14:"R",15:"S",16:"T",17:"V",18:"W",19:"Y",20:"*" ,21 :".",22 : "X"}



def seq2vec(sequence):
    seq_vec = np.zeros((len(sequence), 23))
    i=0
    for ch in sequence:

        seq_vec[i,AA2NUM[ch]]=1
        i+=1
    return seq_vec
def update_concensus(current, seed):
    #print current.shape, seed.shape
    return np.add(current,seed)
def vec2seq(vec):
    maxs = np.argmax(vec, axis=1)
    con_seq = ""
    for x in maxs:
        con_seq = con_seq + NUM2AA[x]
    return con_seq

def concensus_machine_CM(sequences):
    concensus=np.zeros((len(sequences[0]),23))
    for sequence in sequences:
        concensus=update_concensus(concensus,seq2vec(sequence))

    return vec2seq(concensus)


def get_con_seq(sequences):

    if len(sequences)==1:
        return sequences[0]['CDR3']
    concensus = np.zeros((len(sequences[0]['CDR3']), 23))

    for sequence in sequences:
        #print len(sequence['CDR3']), sequence['CDR3']
        concensus = update_concensus(concensus, seq2vec(sequence["CDR3"]))
    return vec2seq(concensus)
if __name__ == '__main__':
    print (concensus_machine_CM(["AAGDDFWSGYSV",
                                "ARGYDFWSGYCY",
                                "ARGYDFWSGYSY",
                                "ARGYDFWSGYQN",
                                "AAGYDFWSGYYF",
                                "AGDDFWSGYFGF"]))




    '''
    [[ 4.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  3.  0.  0.
   0.  0.]
 [ 0.  1.  2.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]]
    '''

