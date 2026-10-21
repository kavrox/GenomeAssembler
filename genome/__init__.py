"""In-Memory Genomic Sequence Assembler & Mutation Analyzer.

The package is split into four parts, one per team member:

    Part 1  sequence, io, kmers, report             String ADT, input/output, k-mers, reports
    Part 2  graph                                   Graph ADT / De Bruijn graph + cleaning
    Part 3  euler, assembler, simulate, benchmark   Eulerian traversal and assembly pipeline
    Part 4  suffix_tree, mutations                  Generalized suffix tree + mutation analysis

The two halves (assembly and mutation analysis) only exchange plain strings.
"""

__version__ = "1.0.0"
