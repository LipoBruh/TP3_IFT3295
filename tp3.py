from boyer_moore import Alignment  # Import specific items

class Query:
    def __init__(self,query_string,seed_size):
       self.query_string = query_string
       self.seed_size = seed_size
       self.kmer_list = [query_string[i:i+seed_size] for i in range(len(query_string) - seed_size + 1)]

    def nth_kmer(self,n):
        if (n < len(self.kmer_list)):
            return self.kmer_list[n]
        return None

    def reverse_complement(self):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(self.query_string))

    def find_seeds(database):
        pass



class Database:
    def __init__(self):
       pass



    def index(query):
        pass


