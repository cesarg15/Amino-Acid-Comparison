#This will import the parser function of BioPython
from Bio import SeqIO
#Keeps a dictionary of integers to add to
from collections import defaultdict
#Will be used to create a bar graph comparing amino acid content
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

#Insert fasta file here
input_file1 = 'lettuce_genome.fa'
input_file2 = 'aphid_genome.fa'

input_files = [input_file1, input_file2]

def count_amino_acid(input_file):
    total_amino_acid_count = defaultdict(int)

    for seq_record in SeqIO.parse(input_file, 'fasta'):
        dna_seq = seq_record.seq
        rna_seq = str(dna_seq.transcribe())

        #print(seq_record.id)
        #print('DNA Sequence: ', seq_record.seq)
        #print(len(seq_record))
        #print('-------------')
        #print('RNA Sequence: ', rna_seq)
        #print(len(rna_seq))
        #print('-------------')

        #Dictionary for codons
        codons = {'UUU':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu','CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu','AUU':'Ile','AUC':'Ile','AUA':'Ile','AUG':'Met','GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val','UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser','AGU':'Ser','AGC':'Ser','CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro','ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala','UAU':'Tyr','UAC':'Tyr','UAA':'Stop','UAG':'Stop','UGA':'Stop','CAU':'His','CAC':'His','CAA':'Gln','CAG':'Gln','AAU':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAU':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu','UGU':'Cys','UGC':'Cys','UGG':'Trp','CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg','AGA':'Arg','AGG':'Arg','GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}

        amino_acid_counts = defaultdict(int)

        codon_list = []
        for codon_start_index in range(0,len(rna_seq),3):
            codon_list.append(rna_seq[codon_start_index:codon_start_index+3])

        for codon,amino_acid in codons.items():
            amino_acid_counts[amino_acid] += codon_list.count(codon)

        #Combines total values for each amino_acid and adds them to total amino acid count
        for amino_acid,count in amino_acid_counts.items():
            total_amino_acid_count[amino_acid] += count

        #print(amino_acid_counts)
        #print('**************')
    result = (dict(total_amino_acid_count))
    return result

input_file1_dict = count_amino_acid(input_file1)
input_file2_dict = count_amino_acid(input_file2)
print('Input 1 Amino Acids :', input_file1_dict)
print('-------------')
print('Input 2 Amino Acids :', input_file2_dict)

#Will convert the dictionaries into arrays of values
input_file1_list = list(input_file1_dict.values())
input_file1_array = np.array(input_file1_list)
input_file2_list = list(input_file2_dict.values())
input_file2_array = np.array(input_file2_list)

print('**************')
paired_ttest = stats.ttest_rel(input_file1_array, input_file2_array)
print('Paired t-test results :', paired_ttest)

X = np.arange(len(input_file1_dict))
ax = plt.subplot(111)
ax.bar(X, input_file1_dict.values(), width=0.2, color='b', align='center')
ax.bar(X+0.2, input_file2_dict.values(), width=0.2, color='g', align='center')
ax.legend(('Input 1','Input 2'))
plt.xticks(X, input_file1_dict.keys())
plt.title("Amino Acid Homogeneity", fontsize=17)
plt.show()
