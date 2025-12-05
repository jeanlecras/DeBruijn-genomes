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
    
    
    def compare(self, fasta_file):
        headers = []
        color_counts = []
        nb_kmer_per_genome = []
        with open(fasta_file, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    prev_line_end = ""
                    header = line[1:].split(" ")[0]
                    headers.append(header)
                    color_counts.append([0]*self.color_count)
                    nb_kmer_per_genome.append(0)
                else:
                    extended_line = prev_line_end+line
                    nb_kmer = len(extended_line)-self.k+1
                    nb_kmer_per_genome[-1] += nb_kmer                
                    for pos in range(nb_kmer):
                        kmer = extended_line[pos:pos+self.k]
                        if kmer in self.table:                           
                            kmer_colors = self.table[kmer]
                            for color in kmer_colors:
                                color_counts[-1][color-1] += 1
                    prev_line_end = line[-self.k+1:]
        
        output_lines = []
        for genome in range(len(nb_kmer_per_genome)):
            scores = [f"{color_count / nb_kmer_per_genome[genome]:.4f}"for color_count in color_counts[genome]]
            line = headers[genome]+"\t" + "\t".join(scores)
            output_lines.append(line)
    
        return "\n".join(output_lines)