
class Alignment:
    def __init__(self):
        self.pattern = ""
        self.algo = BoyerMoore()
        self.aligments = {}

    def set_alignment_algorithm(self,algo):
        self.algo = algo

    def set_pattern(self,pattern):
        self.pattern = pattern

    def align_one(self,pattern,text):
        answers = self.algo.search(text,pattern)
        self.aligments[text] = answers

    def align_all(self,pattern, texts):
        for text in texts:
            self.align_one(pattern,text)






#Source of the python pseudocode :
#Reworked for personnal clarity 
#https://medium.com/@siddharth.21/the-boyer-moore-string-search-algorithm-674906cab162
#
class BoyerMoore:
    def __init__(self,pattern=""):
       self.pattern=pattern
       self.bad_char_table = self.preprocess_bad_character_rule(self.pattern)
       self.good_suffix_table = self.preprocess_good_suffix_rule(self.pattern) 
       pass


    def preprocess_bad_character_rule(self,pattern):
        bad_char_table = {}
        #
        for i in range(len(pattern)):
            bad_char_table[pattern[i]] = len(pattern) - i - 1
        #
        self.bad_char_table = bad_char_table



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





    def search(self,text,pattern):
        answers = []
        #
        if (pattern == self.pattern): #we recycle the tables if same pattern
            self.pattern = pattern
            self.preprocess_bad_character_rule(pattern) #shift table
            self.preprocess_good_suffix_rule(pattern) #
        #
        index_text = len(pattern) - 1
        #
        while (index_text<len(text)-1):

            offset_pattern = 0
            index_pattern = len(pattern)-offset_pattern-1
            #
            while ( (offset_pattern < len(pattern)) and (pattern[index_pattern] == text[index_text-offset_pattern]) ):
                offset_pattern+=1
            #
            if offset_pattern == len(pattern):
                answers.append(index_text+offset_pattern+1) #ajouter une solution 
                index_text+=1
                continue
            else:
                skip1 = self.bad_char_rule(text,index_text,offset_pattern)
                skip2 = self.good_suffix_rule(offset_pattern)
                skip = max(skip1,skip2)
                index_text+= skip
        #
        return answers


    
    #logique découplée pour clarté
    #retourne le shift de l'index du texte
    def bad_char_rule(self,text,index_text,offset_pattern):
        shift = self.bad_char_table.get(text[index_text+offset_pattern],len(self.pattern)) #bad char rule
        if shift == 0 : shift = len(self.pattern)-1
        skips = shift-offset_pattern
        return skips



    def good_suffix_rule(self,offset_pattern):
        skips = self.good_suffix_table[offset_pattern]
        return skips
