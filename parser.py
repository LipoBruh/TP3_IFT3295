class Parser:
    def __init__(self, path = "tRNAs.fasta"):
       self.path = path
       self.lines = []
       self.readFile()

    def readFile(self):
        path = self.path 
        with open(path, 'r') as file:
            lines = [line.strip() for line in file]  # Read lines into a list, stripping newline characters
            self.lines = lines
            return lines

    def get_nth_line(self,n):
        if n<len(self.lines):
            return self.lines[n]
        return None


#decoupled logic for clarity

class Database:
    def __init__(self, path="./tRNAs.fasta"):
        self.parser = Parser(path)
        self.reads=[]
        self.size=0
        self.get_all_reads()

    def get_nth_read(self,n):
        counter=0
        index=0
        symbol = ">"
        #
        for line in self.parser.lines:
            if symbol in line:
                counter+=1
            if counter==n:
                return self.parser.lines[index+1] #returns the rna seq following the nth ">"
            index+=1
        return None
    
    def get_all_reads(self):
        self.size=0
        list=[]
        index=0
        symbol = ">"
        #
        for line in self.parser.lines:
            read = Read()
            if symbol in line:
                read.set_header(line)
                read.set_seq(self.parser.lines[index+1]) #returns the rna seq following the nth ">"
                self.size+=len(self.parser.lines[index+1])
                list.append(read)
            index+=1
        self.reads = list


class Read:
    def __init__(self,path ="tRNAs.fasta",sequence="",header=">"):
        self.path = path
        self.header = header
        self.sequence = sequence

    def set_seq(self,seq):
        self.sequence = seq

    def set_header(self,header):
        self.header = header

    
