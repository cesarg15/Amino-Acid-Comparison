"I was unsure how to create actual test code for my project so instead I opted on making a literal test of the core process of it"
"Here I've removed the requirement for a fasta file to allow more control over the input"
"Although fasta files as input is an important part of my code, the BioPython.parse function is proven to work so I'm not too worried about testing it"
"The input is now a simple string of RNA nucleotides with a set number of codons that should be counted"

#Allows use of a default dictionary of integers to add to
from collections import defaultdict

rna_seq = 'AUGUUUAAUAAAUGCUGGUAA'

#Dictionary for codons
codons = {'UUU':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu','CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu','AUU':'Ile','AUC':'Ile','AUA':'Ile','AUG':'Met','GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val','UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser','AGU':'Ser','AGC':'Ser','CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro','ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala','UAU':'Tyr','UAC':'Tyr','UAA':'Stop','UAG':'Stop','UGA':'Stop','CAU':'His','CAC':'His','CAA':'Gln','CAG':'Gln','AAU':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAU':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu','UGU':'Cys','UGC':'Cys','UGG':'Trp','CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg','AGA':'Arg','AGG':'Arg','GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}
#Holds the amino acid count
amino_acid_counts = defaultdict(int)
#Holds a list if indexed codons
codon_list = []
#Indexes the codons in RNA sequence
for codon_start_index in range(0,len(rna_seq),3):
    codon_list.append(rna_seq[codon_start_index:codon_start_index+3])
#Counts the number of each amino acid from the indexed codon list
for codon,amino_acid in codons.items():
    amino_acid_counts[amino_acid] += codon_list.count(codon)

print(amino_acid_counts)
"Result should be: {'Phe': 1, 'Met': 1, 'Stop': 1, 'Asn': 1, 'Lys': 1, 'Cys': 1, 'Trp': 1}, with the other amino acids being valued at 0"
"The input RNA sequence can be changed for further testing"
"Regardless, this demonstrates that my code is splitting up nucleotides in sets of three and counting only those"
"This excludes any potential overlap between codons being read"
"The rest of my code only adds these dictionaries together and then either graphs that data or converts it into another type"
