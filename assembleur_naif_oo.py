import pickle #pour enregistrer les objets en binaire

class Graph:
    
    def __init__(self, cdbg_file=None, k=17):
        if cdbg_file != None:
            with open(cdbg_file, 'rb') as file:
                loaded_obj = pickle.load(file)
                self.__dict__.update(loaded_obj.__dict__)
        else:
            self.table = {}
            self.k = k
            self.color_count = 0
        
    def __str__(self):
        nb_edges = sum(len(node["prev_kmers"])+len(node["next_kmers"]) for node in self.table.values())
        return f"{len(self.table)} {self.k}mers, {nb_edges} edges, {self.color_count} genomes"
    
    # un génome par fichier
    def add_genome(self, fasta_file):
        self.color_count += 1
        with open(fasta_file,"r") as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):
                    # premier noeud
                    kmer = line[:self.k]
                    next_kmer = line[1:self.k+1]
                    if kmer in self.table:
                        self.table[kmer]["next_kmers"].add(next_kmer)
                        self.table[kmer]["colors"].add(self.color_count)
                    else:
                        self.table[kmer] = {"prev_kmers":set(), "next_kmers":{next_kmer}, "colors":{self.color_count}}
                    
                    for i in range(1,len(line)-self.k):
                        prev_kmer = line[i-1:i-1+self.k]
                        kmer = line[i:i+self.k]
                        next_kmer = line[i+1:i+1+self.k]
                        if kmer in self.table:
                            self.table[kmer]["prev_kmers"].add(prev_kmer)
                            self.table[kmer]["next_kmers"].add(next_kmer)
                            self.table[kmer]["colors"].add(self.color_count)
                        else:
                            self.table[kmer] = {"prev_kmers":{prev_kmer}, "next_kmers":{next_kmer}, "colors":{self.color_count}}
                    
                    #dernier noeud
                    prev_kmer = line[len(line)-self.k-1:len(line)-1]
                    kmer = line[len(line)-self.k:len(line)]
                    if kmer in self.table:
                        self.table[kmer]["prev_kmers"].add(prev_kmer)
                        self.table[kmer]["colors"].add(self.color_count)
                    else:
                        self.table[kmer] = {"prev_kmers":{prev_kmer}, "next_kmers":set(), "colors":{self.color_count}}
                    
    
    def compare(self, fasta_file):
        colors_counter = {}
        kmer_set = set()
        with open(fasta_file,"r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    print(line) #header
                else:
                    for i in range(len(line)-self.k+1):
                        kmer = line[i:i+self.k]
                        if kmer in self.table:
                            for color in self.table[kmer]["colors"]:
                                if color in colors_counter:
                                    colors_counter[color] +=1
                                else:
                                    colors_counter[color] = 1
                        kmer_set.add(kmer)
        
        for color in colors_counter:
            ratio = colors_counter[color]/len(kmer_set)
            print(f"{color} {ratio}", end="\t")
        print()
        
    def write(self, filepath):
        with open(filepath, 'wb') as file:
            pickle.dump(self, file)
                        
    
graph = Graph()
graph.add_genome("covid_ref_genome.fasta")
print(graph)
graph.write("covid.cdbg")