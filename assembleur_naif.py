import pickle #pour enregistrer les objets en binaire

graph = {}

# un génome par fichier
def add_genome(graph, fasta_file, k):
    graph = {}
    color = 1
    with open(fasta_file,"r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                # premier noeud
                kmer = line[:k]
                next_kmer = line[1:k+1]
                if kmer in graph:
                    graph[kmer]["next_kmers"].add(next_kmer)
                    graph[kmer]["colors"].add(color)
                else:
                    graph[kmer] = {"prev_kmers":set(), "next_kmers":{next_kmer}, "colors":{color}}
                
                for i in range(1,len(line)-k):
                    prev_kmer = line[i-1:i-1+k]
                    kmer = line[i:i+k]
                    next_kmer = line[i+1:i+1+k]
                    if kmer in graph:
                        graph[kmer]["prev_kmers"].add(prev_kmer)
                        graph[kmer]["next_kmers"].add(next_kmer)
                        graph[kmer]["colors"].add(color)
                    else:
                        graph[kmer] = {"prev_kmers":{prev_kmer}, "next_kmers":{next_kmer}, "colors":{color}}
                
                #dernier noeud
                prev_kmer = line[len(line)-k-1:len(line)-1]
                kmer = line[len(line)-k:len(line)]
                if kmer in graph:
                    graph[kmer]["prev_kmers"].add(prev_kmer)
                    graph[kmer]["colors"].add(color)
                else:
                    graph[kmer] = {"prev_kmers":{prev_kmer}, "next_kmers":set(), "colors":{color}}
                
    return graph

def compare(fasta_file, graph, k):
    colors_counter = {}
    kmer_set = set()
    with open(fasta_file,"r") as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                print(line) #header
            else:
                for i in range(len(line)-k+1):
                    kmer = line[i:i+k]
                    if kmer in graph:
                        for color in graph[kmer]["colors"]:
                            if color in colors_counter:
                                colors_counter[color] +=1
                            else:
                                colors_counter[color] = 1
                    kmer_set.add(kmer)
    
    for color in colors_counter:
        ratio = colors_counter[color]/len(kmer_set)
        print(f"{color} {ratio}", end="\t")
    print()
    
def save_graph(graph, filepath):
    with open(filepath, 'wb') as file:
        pickle.dump(graph, file)
        
def load_graph(filepath):
    with open(filepath, 'rb') as file:
        return pickle.load(file)
                    

unitig = add_genome(graph, "covid_ref_genome.fasta", 17)