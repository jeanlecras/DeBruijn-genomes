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
            self.color_count = 0 # le premier génome ajouté a la couleur 1
        
    def __str__(self):
        return f"cDBG of {len(self.table)} {self.k}mers and {self.color_count} genomes"
    
    # un génome par fichier
    def add_genome(self, fasta_file):
        self.color_count += 1
        with open(fasta_file,"r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    pass
                else:
                    for pos in range(len(line)-self.k+1):
                        kmer = line[pos:pos+self.k]                        
                        if kmer in self.table:
                            self.table[kmer].add(self.color_count)
                        else:
                            self.table[kmer] = {self.color_count}
                    
    
    def compare(self, fasta_file):
        """
        Questions :
            Peut on ajouter des génomes au fur et à mesure ?
            Faut il ajouter tous ses génomes et après compresser le graphe, ou ajouter un premier génome de façon minamale, puis un autre génome qui garde une structure minamale ?
            Est-ce la que la recherche de génome dans le graphe doit expliquer la fait que le graphe soit compressé ou alors on a le droit de décompresser le graphe pour faire la recherche ?
        """
        headers = []
        color_counts = []
        nb_kmer_per_genome = []
        with open(fasta_file, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:].split(" ")[0]
                    headers.append(header)
                    color_counts.append([0]*self.color_count)
                    nb_kmer_per_genome.append(0)
                else:
                    nb_kmer = len(line)-self.k+1
                    nb_kmer_per_genome[-1] += nb_kmer
                    for pos in range(nb_kmer):
                        kmer = line[pos:pos+self.k] 
                        if kmer in self.table:
                            kmer_colors = self.table[kmer]
                            for color in kmer_colors:
                                color_counts[-1][color-1] += 1
        
        for genome in range(len(nb_kmer_per_genome)):
            print(headers[genome], *(round(color_count/nb_kmer_per_genome[genome],4) for color_count in color_counts[genome]), sep="\t")
        
    def compress(self):
        """
        crée une nouvelle table destinée à stocker des unitigs
        tant que la table otiginale n'est pas vide:
            récupérer le premier kmer de la table originale
            si le kmer actuel a exactement une séquence qui le précède dans la nouvelle table, et cette séquence a extactement 1 kmer qui la suit dans la table originale (alors il s'agit forcément du kmer actuel):
                ajouter une nouvelle entrée à la nouvelle table dont la clé correspond à la séquence (unitig) suivi du dernier nucléotide de notre kmer et la valeur est celle de la séquence (garde les mêmes couleurs)
                supprimer la séquence de la nouvelle table
            sinon:
                ajouter le kmer à la nouvelle table avec la valeur qu'il a dans la table originale'
            supprimer le kmer de la table originale
        remplacer la table originale par la nouvelle table
                                      
        """
        new_table = {}
        while len(self.table) != 0:
            kmer = next(iter(self.table))
            prev_nodes_kmer = set(unitig for unitig in new_table if unitig.endswith(kmer[:-1]))
            unitig = next(iter(prev_nodes_kmer))
            next_nodes_unitig = set(unitig for unitig in self.table if unitig.startswith(kmer[1:]))
            if len(prev_nodes_kmer) == len(next_nodes_unitig) == 1:
                new_table[unitig+kmer[-1]] = new_table[unitig]
                del new_table[unitig]
            else:
                new_table[kmer] = self.table[kmer]
            del self.table[kmer]
        self.table = new_table
            
        
        
    def write(self, filepath):
        with open(filepath, 'wb') as file:
            pickle.dump(self, file)
                        
    
graph = Graph()
graph.add_genome("B253_GCA_000190495.1.fasta")
graph.add_genome("ECOR37_GCA_002190835.1.fasta")
graph.add_genome("H299_GCA_000176695.2.fasta")
graph.add_genome("IAI39_GCA_000026345.1.fasta")
graph.add_genome("SAKAI_GCA_003028755.1.fasta")
print(graph)
graph.compare("query.fa")
