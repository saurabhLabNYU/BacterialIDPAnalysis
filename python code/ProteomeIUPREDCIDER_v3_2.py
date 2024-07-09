
#BacDisDB engine
# Authors: Saumya Saurabh, Ph.D., Stanford University
# Sakshi Kumar, Santa Clara High School & UC Berkeley
# File data names to be entered in line # 28
# from Bio import SeqIO
# for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#      print(seq_record.id)
#      print(repr(seq_record.seq))
#      print(len(seq_record))
# v2 adds sequence cleaner



import os
import csv
from Bio import SeqIO
import sys, re
import pandas as pd
import numpy as np
import iupred2a.iupred2a_lib as iupred2a_lib
import numpy as np
from itertools import groupby
from operator import itemgetter
import localcider
from localcider.sequenceParameters import SequenceParameters

filepath = "C:/Users/sslab_power/Documents/iupred fasta files/ESKAPE/K_UP000008890_2024_05_05.fasta"

def get_longest_consecutive_numbers(numbers):
    idx = max(
        (
            list(map(itemgetter(0), g))
            for i, g in groupby(enumerate(np.diff(numbers)==1), itemgetter(1))
            if i
        ),
        key=len
    )
    return idx[-1]+1-idx[0], idx[0], idx[-1]+1
#returns disordered region length, start residue, end residue
#



records0 = list(SeqIO.parse(filepath, "fasta"))
N1 = int(len(records0))
#adding sequence cleaner here
right_sequences = []
for i in range(0,N1):
    #print(i)
    S = str(records0[i].seq)
    not_AA = ["X","U"]
    if any(x in S for x in not_AA):
        print("Found!")
        print(records0[i])
    else:
        #print("Not found!")
        right_sequences.append(records0[i])
        SeqIO.write(right_sequences, filepath, "fasta")

print("sequence cleaner run")
#end sequence cleaner

#This block uses the cleaned up fasta file and defines outputs

records = list(SeqIO.parse(filepath, "fasta"))
N = int(len(records0))
out_list = [];
cols1 = ['num', 'OS', 'iupred', 'res', 'startend', 'dl', 'length']
cols2 = ['num', 'OS','kappa', 'omega', 'Fdp', 'f+', 'f-', 'Fexp','delta','mw','phosites']
N = int(len(records))
dat1 = pd.DataFrame(columns=cols1)
dat2 = pd.DataFrame(columns=cols2)

print("data files initalized")

# longest sequence
for i in range(0,N):
    #print(i)
    S = str(records[i].seq)
    #print(S)
    iupred_disorder_score = (iupred2a_lib.iupred(S))
    # print(iupred_disorder_score[0])
    res = [idx for idx, val in enumerate(iupred_disorder_score[0]
                                         ) if val > 0.5]
    if len(res) > 10:
        dislen = get_longest_consecutive_numbers(res)[0]
        startend = get_longest_consecutive_numbers(res)
    else:
        dislen = 0
        startend = 0
    protlen = len(iupred_disorder_score)
    orgdes = (records[i].description)
    # print(i)
    # print(res)
    #print(orgdes)
    orgstart = orgdes.find('|') + 1
    orgend = orgdes.find('OX', orgstart)
    OS = orgdes[orgstart:orgend]
    orgName = OS
    nam = records[i].name
    object = SequenceParameters(S)
    k = object.get_kappa()
    om = object.get_Omega()
    dp = object.get_fraction_disorder_promoting()
    fp = object.get_fraction_positive()
    fn = object.get_fraction_negative()
    fexp = object.get_fraction_expanding()
    delta = object.get_delta()
    mw = object.get_molecular_weight()
    phosites = object.get_all_phosphorylatable_sites()

    dat1 = dat1.append({'num': i+1, 'OS': OS, 'iupred': iupred_disorder_score, 'res': res, 'startend': startend, 'dl': dislen, 'length': protlen}, ignore_index=True)
    dat2 = dat2.append({'num': i+1, 'OS': OS, 'kappa': k, 'omega': om, 'Fdp': dp, 'f+': fp, 'f-': fn, 'Fexp': fexp, 'delta': delta, 'mw': mw,'phosites': phosites}, ignore_index=True)
   # print(dat)
    #dat1.to_csv('ProteomeCauloCiderProps20210735.csv')
    with pd.ExcelWriter(filepath+'.xlsx') as writer:
        dat1.to_excel(writer, sheet_name='dat1')
        dat2.to_excel(writer, sheet_name='dat2')

print("Analysis done")

