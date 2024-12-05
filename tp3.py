from boyer_moore import Alignments
from parser import Database , Read
import numpy as np
import os

class Pattern:
    def __init__(self, pattern_string, seed_size):
       self.pattern_string = pattern_string
       self.seed_size = seed_size
       self.kmer_list = self.kmer_list()

    def kmer_list(self):
        list = []
        for i in range (0, len(self.pattern_string), self.seed_size):
            if (len(self.pattern_string[i:])>self.seed_size):
                list.append( self.pattern_string[ i : i+self.seed_size ] )
            else:
                list.append(self.pattern_string[i:])
        return list

    def nth_kmer(self,n):
        if (n < len(self.kmer_list)):
            return self.kmer_list[n]
        return None

    def reverse_complement(self):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(self.pattern_string))
    
    def previous_kmer(self,kmer):
        index = -1
        if kmer in self.kmer_list:
            index = self.kmer_list.index(kmer)
        if index>0:
            return self.kmer_list[index-1]
        return None
    
    def next_kmer(self,kmer):
        index = -1
        if kmer in self.kmer_list:
            index = self.kmer_list.index(kmer)
        if index<(len(self.kmer_list)-1):
            return self.kmer_list[index+1]
        return None
    


class PLAST:

    def __init__(self,pattern="ACATCCTTAGCTCAGTAGGATAGAGCAACAGCCTTCTAAGCTGGTGGTCACAGGTTCAAATCCTGTAGGATGTA", seed="11100001111",match_penality=5,mismatch_penality=-5, min_E = 1e3, ss = 4, path = "./tRNAs.fasta" ):
        self.pattern = Pattern(pattern, len(seed))
        self.seed = seed
        self.db = Database(path)
        self.HSP = {}
        self.extended_HSP=[]
        self.match_penality = match_penality
        self.mismatch_penality = mismatch_penality
        self.min_E = min_E
        self.ss = ss


    def run(self,path="new_file.txt"):
            print("_____LOADING_____")
            #
            self.find_all_HSP()
            #
            sorted_list = self.extend_all()
            #
            if len(sorted_list)>0:
                #print(sorted_list)
                fusion = sorted_list[0][1]
                text = fusion.get_metrics()
                #
                output_directory = "output"
                file_path = os.path.join(output_directory, path)
                #
                os.makedirs(output_directory, exist_ok=True)
                with open(file_path, "w") as file:
                    file.write(text)
                #
            else:
                print("No alignment, try with a more flexible seed")




    ##Boyer Moore Exact match algorithm
    def find_all_HSP(self):
        nb_results = 0
        #
        for kmer in self.pattern.kmer_list:
            #
            if kmer not in self.HSP:
                self.HSP[kmer]=Alignments() #initialize
                self.HSP[kmer].set_pattern_seed(kmer,self.seed) 
            #
            for read in self.db.reads:
                self.HSP[kmer].align(read)  #align a kmer with a fasta seq
            nb_results+=self.HSP[kmer].get_size()

        print(str(nb_results) +" HSP found")

#self.HSP[kmer]= alignments
#aligments[read] = [alignment]
#alignment.score -> score

    def get_nth_HSP(self,n):
        liste = list(self.HSP.keys())
        #
        if len(liste)<=n:
            return None
        else:
            key = liste[n]
            return self.HSP[key]

    def show_max(self,n):
        print(self.HSP)
        #
        if len(self.HSP) == 0 : 
            print("There is no Hot Spot")
            return None
        max = self.get_nth_HSP(n).get_max_score()
        if max is None : 
            print("No Alignment For that Kmer")
            return None
        return max

    def extend_max_n(self,n,match_penalty=5,mismatch_penalty=-4,ss=4):
        align = self.show_max(n)
        return self.extend(align,match_penalty,mismatch_penalty,ss)

    def extend(self,align,match_penalty=5,mismatch_penalty=-4,ss=4):
        #
        if align is None: 
            print("No extension possible")
            return None
        text = align.read.sequence
        score = align.score
        kmer = align.pattern
        full_pattern = self.pattern.pattern_string
        #
        index_text_left = align.index
        index_text_right = align.index+len(kmer)
        index_pattern_left = full_pattern.find(kmer)
        index_pattern_right = index_pattern_left+len(kmer)
        #
        truth_left = (index_text_left>0 and index_pattern_left>0)
        truth_right = (index_text_right<len(text)-1 and index_pattern_right<len(full_pattern)-1)
        flag_left = True
        flag_right = True
        #
        #print("Text of size : "+str(len(text)) + " with starting slice ("+str(index_text_left)+","+str(index_text_right)+")")
        #print("Pattern of size : "+str(len(full_pattern)) + " with starting slice ("+str(index_pattern_left)+","+str(index_pattern_right)+")")
        #
        while((truth_right or truth_left) and score>=ss):
            #chars
            char_text_left = text[index_text_left-1]
            char_pattern_left = text[index_pattern_left-1]
            char_text_right = text[index_text_right+1]
            char_pattern_right = text[index_pattern_right+1]
            #print(str(score)+" Left :("+char_text_left+","+char_pattern_left+") right :("+char_text_right+","+char_pattern_right+")")
            #print("Truths: " + str(truth_left)+ " "+ str(truth_right))
            #print("Flags: "+str(flag_left)+" "+str(flag_right)+"      Indexes:" )
            #cases
            if (char_pattern_left==char_text_left and truth_left):
                score+=match_penalty
                index_pattern_left-=1
                index_text_left-=1
            elif(char_pattern_right==char_text_right and truth_right):
                score+=match_penalty
                index_pattern_right+=1
                index_text_right+=1
            elif (truth_left):
                if (score+mismatch_penalty)<ss: 
                    flag_left = False
                else:
                    score+=mismatch_penalty
                    index_pattern_left-=1
                    index_text_left-=1
            elif(truth_right):
                if (score+mismatch_penalty)<ss: 
                    flag_right = False
                else:
                    score+=mismatch_penalty
                    index_pattern_right+=1
                    index_text_right+=1
            else:
                break
            #
            truth_left = (index_text_left>0 and index_pattern_left>0 and flag_left)
            truth_right = (index_text_right<len(text)-1 and index_pattern_right<len(full_pattern)-1 and flag_right)
            #
        print("Ext. Score      : "+str(score))
        print("Ext. Pattern    : "+full_pattern[index_pattern_left:index_pattern_right])
        print("Extended HSP    : "+text[index_text_left:index_text_right])
        return Extended_HSP(align, index_text_left, index_text_right,index_pattern_left,index_pattern_right,self.pattern)
    

    #Used to find kmers that have hot spots in the same sequence
    #We will use the list of KMERs to extend those with the best score
    #  It returns a dict of the form: 
    # {Read : [align1, align2 ...]}
    # Such that for a specific Read / Sequence, we could extract the 2 alignments with the best score to extend them
    def find_common_alignments(self):
        to_extend = {}
        #toutes les paires sans répétition de paires
        for i  in range(0,len(self.pattern.kmer_list),1):
            for j in range(i+1,len(self.pattern.kmer_list),1):
                kmer1 = self.pattern.kmer_list[i]
                kmer2 = self.pattern.kmer_list[j]
                #
                headers1 = self.HSP[kmer1].get_headers() #get all alignment headers for that kmer
                headers2 = self.HSP[kmer2].get_headers()
                #
                common = self.intersection(list(headers1.keys()),list(headers2.keys()))
                if len(common)==0:
                    continue

                #a single Alignments class instance is specific to a kmer
                # HSP[kmer] returns the Alignments class specific to that kmer
                # Alignments[read] returns a list of Alignment class
                # an Alignment class 
                else:
                    for header in common:
                        read = headers1[header]
                        if read in to_extend:
                            to_extend[read].add(kmer1)
                            to_extend[read].add(kmer2)
                        #we use a set to simplify duplicate entries
                        else:
                            to_extend[read]=set((kmer1,kmer2))

            #print(to_extend)
            return to_extend
            #dictionnary with structure : {kmer1: {kmer2: [read1,read2...]}}


    def intersection(self,a, b):
        c = [value for value in a if value in b]
        return c
        
    def kmer_max_score_read(self,kmer,read):
        alignments = self.HSP[kmer]
        alignment = alignments.get_max_read_score(read)
        return alignment
    
    def extend_all(self):
        fusion_dic = {}
        dic = self.find_common_alignments()
        print(str(len(dic))+" reads with 2+ HSP")
        for read in dic:
            fusion_dic[read] = self.extend_read_kmer(read,dic[read]) #should return a Fusion_HSP object
        #
        sorted_items = sorted(fusion_dic.items(), key=lambda pair: pair[1].e_value(), reverse=False)
        #
        condition = lambda item: item[1].cut_off(self.min_E)
        filtered_list = list(filter(condition,sorted_items))
        #
        return filtered_list # dict of the format {read : Fusion_HSP}


    #will find the best alignment for a kmer in a specific read
    #then will extend each kmers that belongs to that read
    #will then create a fusion if possible and return that fusion
    def extend_read_kmer(self,read,kmer_set):
        extended_list = []
        for kmer in kmer_set:
            alignment = self.kmer_max_score_read(kmer,read)
            extended_hsp = self.extend(alignment,self.match_penality,self.mismatch_penality,self.ss)
            extended_list.append(extended_hsp)
        if len(extended_list)==0:
            print("No extension possible")
            return None
        elif len(extended_list)==1:
            print("Only a single sequence was found and extended, no fusion occurred.")
            return Fusion_HSP(extended_list[0],self.db)
        else:
            fusion = Fusion_HSP(extended_list[0],self.db)
            for extension in extended_list[1:]:
                fusion.fuse_HSP(extension)
            return fusion


"""
"""



class Extended_HSP:
    def __init__(self, alignment, index_text_left,index_text_right, index_pattern_left, index_pattern_right, pattern):
       self.alignment = alignment
       self.index_text_left = index_text_left
       self.index_text_right = index_text_right
       self.index_pattern_left = index_pattern_left
       self.index_pattern_right = index_pattern_right
       self.pattern = pattern

    def get_sequence(self):
       return self.alignment.sequence[self.index_text_left : self.index_text_right]

class Fusion_HSP:
    def __init__(self, initial_HSP,db):
        self.initial_HSP = initial_HSP
        self.read = initial_HSP.alignment.read
        self.pattern = initial_HSP.pattern
        self.db = db
        #
        self.index_text_left = initial_HSP.index_text_left
        self.index_text_right = initial_HSP.index_text_right
        self.index_pattern_left = initial_HSP.index_pattern_left
        self.index_pattern_right = initial_HSP.index_pattern_right
        #
        self.fused_HSP = [initial_HSP]

    def codon_to_amino_acid(self,codon):
        codon_map = {
            # Phe
            "UUU": {"acronym": "Phe", "name": "Phenylalanine", "image": "./images/L-Phenylalanine.png","letter":"F"},
            "UUC": {"acronym": "Phe", "name": "Phenylalanine", "image": "./images/L-Phenylalanine.png","letter":"F"},
            # Leu
            "CUU": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            "CUC": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            "CUA": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            "CUG": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            "UUA": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            "UUG": {"acronym": "Leu", "name": "Leucine", "image": "./images/L-Leucine.png","letter":"L"},
            # Ser
            "AGU": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            "AGC": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            "UCU": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            "UCC": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            "UCA": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            "UCG": {"acronym": "Ser", "name": "Serine", "image": "./images/L-Serine.png","letter":"S"},
            # Tyr
            "UAU": {"acronym": "Tyr", "name": "Tyrosine", "image": "./images/L-Tyrosine.png","letter":"Y"},
            "UAC": {"acronym": "Tyr", "name": "Tyrosine", "image": "./images/L-Tyrosine.png","letter":"Y"},
            # STOP
            "UGA": {"acronym": "STOP", "name": "STOP", "image": "./images/STOP.png","letter":"*"},
            "UAA": {"acronym": "STOP", "name": "STOP", "image": "./images/STOP.png","letter":"*"},
            "UAG": {"acronym": "STOP", "name": "STOP", "image": "./images/STOP.png","letter":"*"},
            # Cys
            "UGU": {"acronym": "Cys", "name": "Cysteine", "image": "./images/L-Cysteine.png","letter":"C"},
            "UGC": {"acronym": "Cys", "name": "Cysteine", "image": "./images/L-Cysteine.png","letter":"C"},
            # Trp
            "UGG": {"acronym": "Trp", "name": "Tryptophan", "image": "./images/L-Tryptophan.png","letter":"W"},
            # Pro
            "CCU": {"acronym": "Pro", "name": "Proline", "image": "./images/L-Proline.png","letter":"P"},
            "CCC": {"acronym": "Pro", "name": "Proline", "image": "./images/L-Proline.png","letter":"P"},
            "CCA": {"acronym": "Pro", "name": "Proline", "image": "./images/L-Proline.png","letter":"P"},
            "CCG": {"acronym": "Pro", "name": "Proline", "image": "./images/L-Proline.png","letter":"P"},
            # His
            "CAU": {"acronym": "His", "name": "Histidine", "image": "./images/L-Histidine.png","letter":"H"},
            "CAC": {"acronym": "His", "name": "Histidine", "image": "./images/L-Histidine.png","letter":"H"},
            # Gln
            "CAA": {"acronym": "Gln", "name": "Glutamine", "image": "./images/L-Glutamine.png","letter":"Q"},
            "CAG": {"acronym": "Gln", "name": "Glutamine", "image": "./images/L-Glutamine.png","letter":"Q"},
            # Arg
            "CGU": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            "CGC": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            "CGA": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            "CGG": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            "AGA": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            "AGG": {"acronym": "Arg", "name": "Arginine", "image": "./images/L-Arginine.png","letter":"R"},
            # Ile
            "AUU": {"acronym": "Ile", "name": "Isoleucine", "image": "./images/L-Isoleucine.png","letter":"I"},
            "AUC": {"acronym": "Ile", "name": "Isoleucine", "image": "./images/L-Isoleucine.png","letter":"I"},
            "AUA": {"acronym": "Ile", "name": "Isoleucine", "image": "./images/L-Isoleucine.png","letter":"I"},
            # Met
            "AUG": {"acronym": "Met", "name": "Methionine", "image": "./images/L-Methionine.png","letter":"M"},
            # Thr
            "ACU": {"acronym": "Thr", "name": "Threonine", "image": "./images/L-Threonine.png","letter":"T"},
            "ACC": {"acronym": "Thr", "name": "Threonine", "image": "./images/L-Threonine.png","letter":"T"},
            "ACA": {"acronym": "Thr", "name": "Threonine", "image": "./images/L-Threonine.png","letter":"T"},
            "ACG": {"acronym": "Thr", "name": "Threonine", "image": "./images/L-Threonine.png","letter":"T"},
            # Asn
            "AAU": {"acronym": "Asn", "name": "Asparagine", "image": "./images/L-Asparagine.png","letter":"N"},
            "AAC": {"acronym": "Asn", "name": "Asparagine", "image": "./images/L-Asparagine.png","letter":"N"},
            # Lys
            "AAA": {"acronym": "Lys", "name": "Lysine", "image": "./images/L-Lysine.png","letter":"K"},
            "AAG": {"acronym": "Lys", "name": "Lysine", "image": "./images/L-Lysine.png","letter":"K"},
            # Val
            "GUU": {"acronym": "Val", "name": "Valine", "image": "./images/L-Valine.png","letter":"V"},
            "GUC": {"acronym": "Val", "name": "Valine", "image": "./images/L-Valine.png","letter":"V"},
            "GUA": {"acronym": "Val", "name": "Valine", "image": "./images/L-Valine.png","letter":"V"},
            "GUG": {"acronym": "Val", "name": "Valine", "image": "./images/L-Valine.png","letter":"V"},
            # Ala
            "GCU": {"acronym": "Ala", "name": "Alanine", "image": "./images/L-Alanine.png","letter":"A"},
            "GCC": {"acronym": "Ala", "name": "Alanine", "image": "./images/L-Alanine.png","letter":"A"},
            "GCA": {"acronym": "Ala", "name": "Alanine", "image": "./images/L-Alanine.png","letter":"A"},
            "GCG": {"acronym": "Ala", "name": "Alanine", "image": "./images/L-Alanine.png","letter":"A"},
            # Asp
            "GAU": {"acronym": "Asp", "name": "Aspartate", "image": "./images/L-Aspartate.png","letter":"D"},
            "GAC": {"acronym": "Asp", "name": "Aspartate", "image": "./images/L-Aspartate.png","letter":"D"},
            # Glu
            "GAA": {"acronym": "Glu", "name": "Glutamate", "image": "./images/L-Glutamate.png","letter":"E"},
            "GAG": {"acronym": "Glu", "name": "Glutamate", "image": "./images/L-Glutamate.png","letter":"E"},
            # Gly
            "GGU": { "acronym": "Gly", "name": "Glycine", "image": "./images/L-Glycine.png","letter":"G" },
            "GGC": { "acronym": "Gly", "name": "Glycine", "image": "./images/L-Glycine.png","letter":"G" },
            "GGA": { "acronym": "Gly", "name": "Glycine", "image": "./images/L-Glycine.png","letter":"G" },
            "GGG": { "acronym": "Gly", "name": "Glycine", "image": "./images/L-Glycine.png","letter":"G" },
        }
        if codon in codon_map:
            return codon_map[codon]["letter"]
        return None


    def fuse_HSP(self,extended):
        #same string?
        if self.read.header != extended.alignment.read.header:
            print("Different Sequences")
            return None
        #intersection?
        # Case 1 :      ____===‾‾‾‾    (L1 < L2 and R1 > L2) or (L2 < L1 and R2 > L1)
        # Case 2 :    _____     _____  L1 R1 < L2    or  L2 R2 < L1   or  
        # Case 3 :     _____=====_____ (L1 R1 < R2  and L1 R1 > L2) or (L2 R2 < R1  and L2 R2 > L1)
        if not ((self.index_text_left < extended.index_text_left and self.index_text_right > extended.index_text_left) 
                or (extended.index_text_left < self.index_text_left and extended.index_text_right > self.index_text_left)):
            print("No Intersection between extensions")
            return None
        else:
            self.fused_HSP.append(extended)
            self.index_text_left = min(self.index_text_left, extended.index_text_left)
            self.index_text_right = max(self.index_text_right, extended.index_text_right)
            self.index_pattern_left = min(self.index_pattern_left, extended.index_pattern_left)
            self.index_pattern_right = max(self.index_pattern_right, extended.index_pattern_right)

    def get_sequence(self):
        text = self.read.sequence[self.index_text_left : self.index_text_right]
        return text
    def get_pattern(self):
        text = self.pattern.pattern_string[self.index_pattern_left : self.index_pattern_right]
        return text


    #works
    def get_score(self,match_penalty=5,mismatch_penalty=-4):
            #
            self.score=0
            sequence = self.get_sequence()
            pattern = self.get_pattern()
            #
            for i in range(0,len(sequence), 1):
                #
                if sequence[i]==pattern[i]:
                    self.score+=match_penalty
                else:
                    self.score+=mismatch_penalty
            return self.score
            

    def bitscore(self):
        self.get_score()
        b = round(   (self.score - np.log(0.176))/np.log(2)   )
        return b
    
    def e_value(self):
        e = self.db.size*len(self.pattern.pattern_string)*(2**(-1*self.bitscore()))
        return e
    
    def cut_off(self,min_E):
        return self.e_value()<min_E
        

    def to_reverse(self):
        dna_to_rna = {
        "A": "U",  # Adenine -> Uracil
        "T": "A",  # Thymine -> Adenine
        "C": "G",  # Cytosine -> Guanine
        "G": "C"   # Guanine -> Cytosine
    }
        return "".join([dna_to_rna[base] for base in reversed(self.get_sequence())])
    

    def to_AA(self):
        rna = self.to_reverse()
        print(rna)
        aa_seq = ""
        for i in range (0,len(rna),3):
            if i+3>len(rna):
                return aa_seq+"..."
            codon = rna[i:i+3]
            print(codon)
            aa = self.codon_to_amino_acid(codon)
            aa_seq+=aa
        return aa_seq




    def get_metrics(self):
        text = (
            "_____RESULTS_____\n"
            f"  Seed used      : {self.initial_HSP.alignment.seed}\n"
            f"  Pattern match  : {self.get_pattern()}\n"
            f"  Sequence match : {self.get_sequence()}\n"
            f"  ID             : {self.read.header}\n"
            f"  Fusion score   : {self.get_score()}\n"
            f"  Bit score      : {self.bitscore()}\n"
            f"  E value        : {self.e_value()}\n"
            f"  Q2.2 AA      (deduced) : {self.to_AA()}\n"
            f"  Q2.2 Codon   (deduced) : {self.to_reverse()}\n"
        )

        print(text)
        return text

#Lazy way to do it
# format :     plast =(self,pattern="...", seed="...",match_penality=5,mismatch_penality=-4, min_E = 0.001, ss = 4, path = "./tRNAs.fasta" )
#seq1
plast = PLAST("AGCGGGGTAGAGGAATTGGTTTACTCATCAGGCTCATGACCTGAAGACTGCAGGTTCGAATCCTGTCCCCGCCT", "11111111111",5,-4, 0.001, 4, 'tRNAs.fasta')
plast.run("q2_seq1")
#seq2
plast = PLAST("ACATCCTTAGCTCAGTAGGATAGAGCAACAGCCTTCTAAGCTGGTGGTCACAGGTTCAAATCCTGTAGGATGTA", "11100000111",5,-4, 0.001, 4, 'tRNAs.fasta')
plast.run("q2_seq2")
#seq3
plast = PLAST("CGCGGAGTAGAGCAGTTTGGTAGCTCGCAAGGCTCATAACCTTGAGGTCACGGGTTCAAATCCTGTCATCCCTA", "11100000111",5,-4, 100, 4, 'tRNAs.fasta')
plast.run("q2_seq3")
#seq4
plast = PLAST("GCATTCTTAGCTCAGCTGGATAGAGCAACAACCTTCTAAGTTGAAGGTCACAGGTTCAAATCCTGTAGGATGCT", "11111111111",5,-4, 0.001, 4, 'tRNAs.fasta')
plast.run("q2_seq4")