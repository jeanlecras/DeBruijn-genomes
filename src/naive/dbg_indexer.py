#!/usr/bin/env python3
import argparse
import time
import sys

import build as build_module
import query as query_module


def cmd_build(args):
    """
    dbg_indexer.py build -i sequences_file_list -k kmer_size -o output_file
    """
    start_build = time.time()

    # Chargement de la liste des génomes
    with open(args.i) as f:
        genome_files = [line.strip() for line in f if line.strip()]

    # Construction du graph
    graph = build_module.Graph(k=args.k)

    for genome in genome_files:
        graph.add_genome(genome)

    end_build = time.time()
    print(f"OUT TIME_BUILD: {end_build - start_build:.4f}")

    # Sérialisation
    start_ser = time.time()
    graph.write(args.o)
    end_ser = time.time()
    print(f"OUT TIME_SERIALISATION: {end_ser - start_ser:.4f}")


def cmd_query(args):
    """
    dbg_indexer.py query -q query_file.fa -i cdbg -o output_file
    """
    start_deser = time.time()

    # chargement du graph
    graph = query_module.Graph(cdbg_file=args.i)

    end_deser = time.time()
    print(f"OUT TIME_DESERIALISATION: {end_deser - start_deser:.4f}")

    start_query = time.time()

    result_str = graph.compare(args.q)

    with open(args.o, "w") as out:
        out.write(result_str)

    end_query = time.time()
    print(f"OUT TIME_QUERY: {end_query - start_query:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="De Bruijn Graph Indexer (build & query)"
    )

    subparsers = parser.add_subparsers(
        title="commands",
        dest="command"
    )

    # build
    build_parser = subparsers.add_parser(
        "build",
        help="Construire un cDBG à partir d'une liste de génomes"
    )
    build_parser.add_argument("-i", required=True,
                              help="Fichier listant les chemins FASTA")
    build_parser.add_argument("-k", type=int, required=True,
                              help="Taille des k-mers")
    build_parser.add_argument("-o", required=True,
                              help="Fichier de sortie du cDBG sérialisé")
    build_parser.set_defaults(func=cmd_build)

    # query 
    query_parser = subparsers.add_parser(
        "query",
        help="Interroger un cDBG avec des séquences FASTA"
    )
    query_parser.add_argument("-q", required=True,
                              help="Fichier FASTA contenant les requêtes")
    query_parser.add_argument("-i", required=True,
                              help="Fichier contenant le cDBG sérialisé")
    query_parser.add_argument("-o", required=True,
                              help="Fichier de sortie")
    query_parser.set_defaults(func=cmd_query)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
