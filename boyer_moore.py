#decoupled logic for clarity
from parser import Read

class Alignments:
    def __init__(self, pattern=""):
        self.pattern = pattern #kmer
        self.algo = BoyerMoore()
        self.alignments = {} #list of alignment class items


    def set_pattern(self,pattern):
        self.seed = "1"*len(pattern)
        self.pattern = pattern
        self.alignments={}
        return True
    
    def set_pattern_seed(self,pattern,seed):
        if len(seed)<len(pattern): return None
        self.seed = seed[:len(pattern)]
        self.pattern = pattern
        self.alignments={}
        return True


    def align(self,read):
        results = self.algo.search(read,self.pattern,self.seed)
        if len(results)>0:
            self.alignments[read] = results

    def get_max_score(self):
        #
        if len(self.alignments)==0:
            return None
        max = None
        #
        for read in self.alignments:
            #
            for align in self.alignments[read]:
                
                if align is None : 
                    continue
                if max is None:
                    max = align
                elif align.score>max.score:
                    max = align
        print("...")
        print("Full read       : "+max.read.sequence)
        print("Score           : "+str(max.score))
        print("Seed            : "+self.seed)
        print("Kmer match      : "+self.pattern)
        print("Sequence match  : "+max.get_match_sequence())
        return max
    
    def get_max_read_score(self,read):
        #
        if len(self.alignments)==0:
            return None
        max = None
        #
        for align in self.alignments[read]:
                if align is None : 
                    continue
                if max is None:
                    max = align
                elif align.score>max.score:
                    max = align
        print("...")
        print("Read            : "+max.read.sequence)
        print("Score           : "+str(max.score))
        print("Seed            : "+self.seed)
        print("Kmer match      : "+self.pattern)
        print("Sequence match  : "+max.get_match_sequence())
        return max

    
    #a single Alignments class instance is specific to a kmer
    # HSP[kmer] returns the Alignments class specific to that kmer
    # Alignments[read] returns a list of Alignment class
    # an Alignment class 
    def get_headers(self):
        headers = {}
        for read in self.alignments:
            headers[read.header] = read
        return headers
    
    def get_size(self):
        size = 0
        for e in self.alignments:
            size+=len(self.alignments[e])
        return size

class Alignment:
        def __init__(self, read,pattern,index,seed):
            self.index = index
            self.pattern = pattern
            self.read = read
            self.seed = seed
            self.score = self.calculate_score(5,-4)

        def calculate_score(self,match_penalty,mismatch_penalty):
            #
            self.score=0
            index_pattern=0
            for i in range(self.index, min(self.index+len(self.pattern),len(self.read.sequence)), 1):
                #
                if self.pattern[index_pattern]==self.read.sequence[i]:
                    self.score+=match_penalty
                else:
                    self.score+=mismatch_penalty
                index_pattern+=1
            return self.score
        
        def get_match_sequence(self):
            return self.read.sequence[self.index : self.index+len(self.pattern)]
            





#Source of the python pseudocode :
#Reworked for personnal clarity 
#https://medium.com/@siddharth.21/the-boyer-moore-string-search-algorithm-674906cab162
#
class BoyerMoore:
    def __init__(self,pattern=""):
       self.pattern=pattern
       self.bad_char_table = self.preprocess_bad_character_rule(self.pattern)
       self.good_suffix_table = self.preprocess_good_suffix_rule(self.pattern) 



    def preprocess_bad_character_rule(self,pattern):
        bad_char_table = {}
        #
        for i in range(len(pattern)):
            bad_char_table[pattern[i]] = len(pattern) - i - 1
        #
        self.bad_char_table = bad_char_table


    #errors in the preprocessing so omitted for now
    #https://www.geeksforgeeks.org/boyer-moore-algorithm-good-suffix-heuristic/ heuristics for the good suffix rule
    def preprocess_good_suffix_rule(self, pattern):
        #border = substring that can be both a prefix / suffix of the pattern
        pattern_length = len(pattern)
        border_positions = [0] * (pattern_length + 1)  # Positions of borders, discarded later 
        shift_distances = [0] * (pattern_length + 1)  # Distances to shift, return value

        # Preprocessing for strong suffix
        current_index, border_index = pattern_length, pattern_length + 1
        border_positions[current_index] = border_index  # Border beyond last character

        while current_index > 0:
            # Update shift and border positions for mismatches
            while border_index <= pattern_length and pattern[current_index - 1] != pattern[border_index - 1]:
                if shift_distances[border_index] == 0:  shift_distances[border_index] = border_index - current_index # Set shift if not already set
                border_index = border_positions[border_index]  # Move to the next widest border

            current_index -= 1
            border_index -= 1
            border_positions[current_index] = border_index  # Store new border position

        # Preprocessing for prefix case
        border_index = border_positions[0]
        for index in range(pattern_length + 1):
            if shift_distances[index] == 0:  # If shift not set, use prefix-based shift
                shift_distances[index] = border_index
            if index == border_index:  # Update border index to next widest
                border_index = border_positions[border_index]

        self.good_suffix_table = shift_distances  # Save the shift distances





    def search(self,read,pattern,seed):
        text = read.sequence
        answers = []
        self.pattern = pattern
        self.preprocess_bad_character_rule(pattern) #shift table
        self.preprocess_good_suffix_rule(pattern) #
        #
        index_text = len(pattern) - 1
        #
        while (index_text<len(text)-1):

            offset_pattern = 0
            #index_pattern = len(pattern)-offset_pattern-1
            #index = index_text-offset_pattern
            #    
            while ( (offset_pattern < len(pattern)) and (pattern[len(pattern)-offset_pattern-1] == text[index_text-offset_pattern]) or seed[len(pattern)-offset_pattern-1] == "0"):
                offset_pattern+=1
                
            #
            if offset_pattern == len(pattern):
                answers.append(Alignment(read, pattern,(index_text+1-len(self.pattern)),seed))#ajouter une solution 
                index_text+=1
                continue
            else:
                skip1 = self.bad_char_rule(offset_pattern)
                #skip2 = self.good_suffix_rule(offset_pattern)
                skip = max(skip1,1)
                index_text+= skip
        return answers


    
    #logique découplée pour clarté
    #retourne le shift de l'index du texte
    def bad_char_rule(self,offset_pattern):
        #
        index_pattern = len(self.pattern)-offset_pattern-1
        missmatch_char = self.pattern[index_pattern]
        #
        skip = self.bad_char_table.get(missmatch_char,len(self.pattern)) #return bad char or len of pattern if not in table
        return skip



    def good_suffix_rule(self,offset_pattern):
        skips = self.good_suffix_table[offset_pattern]
        return skips

