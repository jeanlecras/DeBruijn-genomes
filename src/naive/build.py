import pickle

class Graph:
    
    def __init__(self, cdbg_file:str=None, k:int=31):
        """
        Initialise an empty De Bruijn graph or load an existing De Bruijn graph from a cdbg file

        Parameters
        ----------
        cdbg_file : str, optional
            A cdbg file
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
    
                            
    def add_genome(self, fasta_file: str):
        """
        Add a genome into the De Bruijn graph and marks its kmers with its color

        Parameters
        ----------
        fasta_file : str
            A fasta file containing a single genome with one or more contigs
        """
        self.color_count += 1 # incrémente le conteur de génomes
        with open(fasta_file,"r") as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    prev_line_end = ""
                else:
                    extented_line = prev_line_end + line #prend en compte les kmers à cheval sur 2 lignes 
                    for pos in range(len(extented_line)-self.k+1):
                        kmer = extented_line[pos:pos+self.k]                        
                        if kmer in self.table:
                            self.table[kmer].add(self.color_count)
                        else:
                            self.table[kmer] = {self.color_count} # ces kmers ne se trouvent dans aucun autre génome
                    prev_line_end = line[-self.k+1:] #les derniers nucléotides de la ligne qui serviront  former les kmers de la ligne suivante
                                        
        
        
    def write(self, filepath: str):
        """
        Save the python Graph object into a binary file (cdbg file)

        Parameters
        ----------
        filepath : str
            Name of the cdbg file
        """
        with open(filepath, 'wb') as file:
            pickle.dump(self, file)