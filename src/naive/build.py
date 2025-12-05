import pickle

class Graph:
    
    def __init__(self, cdbg_file=None, k=31):
        if cdbg_file != None:
            with open(cdbg_file, 'rb') as file:
                loaded_obj = pickle.load(file)
                self.__dict__.update(loaded_obj.__dict__)
        else:
            self.table = {}
            self.k = k
            self.color_count = 0 # le premier génome ajouté a la couleur 1
        
    def __str__(self):
        return f"cDBG of {len(self.table)} {self.k}mers and {self.color_count} genomes"
    
    
    # un génome pr fichier                        
    def add_genome(self, fasta_file):
        self.color_count += 1
        with open(fasta_file,"r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    prev_line_end = ""
                else:
                    extented_line = prev_line_end + line
                    for pos in range(len(extented_line)-self.k+1):
                        kmer = extented_line[pos:pos+self.k]                        
                        if kmer in self.table:
                            self.table[kmer].add(self.color_count)
                        else:
                            self.table[kmer] = {self.color_count}
                    prev_line_end = line[-self.k+1:]
                                        
        
        
    def write(self, filepath):
        with open(filepath, 'wb') as file:
            pickle.dump(self, file)