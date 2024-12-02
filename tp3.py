from boyer_moore import Alignments
from parser import Database , Read
import numpy as np

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

    def __init__(self,pattern="ACATCCTTAGCTCAGTAGGATAGAGCAACAGCCTTCTAAGCTGGTGGTCACAGGTTCAAATCCTGTAGGATGTA", seed="11100001111",match_penality=5,mismatch_penality=-5, min_E = 4, ss = 1e-3, path = "./tRNAs.fasta" ):
        self.pattern = Pattern(pattern, len(seed))
        self.seed = seed
        self.db = Database(path)
        self.HSP = {}
        self.extended_HSP=[]
        self.match_penality = match_penality
        self.mismatch_penality = mismatch_penality
        self.min_E = min_E
        self.ss = ss


    def run(self):
            print("_____LOADING_____")
            #
            self.find_all_HSP()
            #
            sorted_list = self.extend_all()
            #
            if len(sorted_list)>0:
                #print(sorted_list)
                fusion = sorted_list[0][1]
                fusion.get_metrics()
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
        text = max.read.sequence
        index = max.index
        return max
        print("Text match : " + text[index:index+11])

    def extend_max_n(self,n,match_penalty=5,mismatch_penalty=-4,E=4):
        align = self.show_max(n)
        return self.extend(align,match_penalty,mismatch_penalty,E)

    def extend(self,align,match_penalty=5,mismatch_penalty=-4,E=4):
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
        while((truth_right or truth_left) and score>=E):
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
                if (score+mismatch_penalty)<E: 
                    flag_left = False
                else:
                    score+=mismatch_penalty
                    index_pattern_left-=1
                    index_text_left-=1
            elif(truth_right):
                if (score+mismatch_penalty)<E: 
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
        return sorted_items # dict of the format {read : Fusion_HSP}


    #will find the best alignment for a kmer in a specific read
    #then will extend each kmers that belongs to that read
    #will then create a fusion if possible and return that fusion
    def extend_read_kmer(self,read,kmer_set):
        extended_list = []
        for kmer in kmer_set:
            alignment = self.kmer_max_score_read(kmer,read)
            extended_hsp = self.extend(alignment,self.match_penality,self.mismatch_penality,self.min_E)
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
        
    def get_metrics(self):
        print("_____RESULTS_____")
        print("  Seed used      : " + self.initial_HSP.alignment.seed)
        print("  Pattern match  : "+ self.get_pattern())
        print("  Sequence match : "+ self.get_sequence())
        print("  ID             : "+self.read.header)
        print("  Fusion score   : "+str(self.get_score()))
        print("  Bit score      : "+str(self.bitscore()))
        print("  E value        : "+str(self.e_value()))

"""
CGAATCCTTCTAGA
AGGTGAGGGTTC
AGGTCTATGGCTTAAGGGTAGAGTGTTGGTTTCCAAAACCAC AGGTGAGGGTTCGAATCCTTCTAGA CCTG

"""

"""
pattern CGTAGTCGGCTAACCAGCATAACGCTTGTAAACGTAAGAGCCC

AGGTGAGGGTTCGAATCCTTCTAGA
"""