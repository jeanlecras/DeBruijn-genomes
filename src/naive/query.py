import pickle

class Graph:
    
    def __init__(self, cdbg_file:str=None, k:int=31):
        """
        Initialise an empty De Bruijn graph or load an existing De Bruijn graph from a cdbg file

        Parameters
        ----------
        cdbg_file : str, optional
            The address of a cdbg file
        k : int, optional
            The size of the kmers. Must be greater than 16. The default is 31
        """
        if k < 17:
            raise ValueError(f"k-mer size must be >= 17 (received: {k})")
            
        if cdbg_file != None:
            with open(cdbg_file, 'rb') as file:
                loaded_obj = pickle.load(file)
                self.__dict__.update(loaded_obj.__dict__)
        else:
            self.table = {} # un dictionnaire dont les clés sont les kmers et les valeurs sont des ensembles de couleurs de génomes
            self.k = k # taille des kmers
            self.color_count = 0 # nombre de génomes, représente la couleur qu'on associe au dernier génome ajouté
        
    def __str__(self):
        return f"cDBG of {len(self.table)} {self.k}mers and {self.color_count} genomes"
    
    
    def compare(self, fasta_file: str) -> str:
        """
        Compare genomes to the genomes of the graph to get for each of them, how much kmers in common do they share with the genomes of the graph

        Parameters
        ----------
        fasta_file : str
            A fasta file containing one or more genomes

        Returns
        -------
        TYPE: str
            The alignment rate of query genomes (rows) on graph genomes (columns) 
        """
        headers = [] #liste contenant le début des headers de chaque génome
        color_counts = [] #Matrice des nombres de couleurs en commun avec en colonnes les couleurs (génomes du graphe) et en ligne les génomes query
        nb_kmer_per_genome = [] #nombre de kmer de chacun des génomes query
        with open(fasta_file, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    prev_line_end = ""
                    header = line[1:].split(" ")[0]
                    headers.append(header)
                    color_counts.append([0]*self.color_count) #les comptes de chaque couleur sont initialisés à 0
                    nb_kmer_per_genome.append(0)
                else:
                    extended_line = prev_line_end+line #prend en compte les kmers à cheval sur 2 lignes 
                    nb_kmer = len(extended_line)-self.k+1
                    nb_kmer_per_genome[-1] += nb_kmer                
                    for pos in range(nb_kmer):
                        kmer = extended_line[pos:pos+self.k]
                        if kmer in self.table:                           
                            kmer_colors = self.table[kmer]
                            for color in kmer_colors:
                                color_counts[-1][color-1] += 1
                    prev_line_end = line[-self.k+1:] #les derniers nucléotides de la ligne qui serviront  former les kmers de la ligne suivante
        
        rows = []
        for genome in range(len(nb_kmer_per_genome)):
            scores = [f"{color_count / nb_kmer_per_genome[genome]:.4f}" for color_count in color_counts[genome]]
            line = headers[genome]+"\t" + "\t".join(scores)
            rows.append(line)
            tab = "\n".join(rows)
    
        return tab