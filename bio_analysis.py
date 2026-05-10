"""
Bioinformatics Analysis Script
Provides tools for DNA sequence analysis, genomic data processing, and basic statistics.
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter


class SequenceAnalyzer:
    """
    A class to analyze DNA sequences and perform bioinformatics calculations.
    """
    
    def __init__(self, sequence):
        """
        Initialize with a DNA sequence.
        
        Args:
            sequence (str): DNA sequence string (ATCG)
        """
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
    
    def get_length(self):
        """Return the length of the sequence."""
        return self.length
    
    def get_gc_content(self):
        """
        Calculate GC content (percentage of G and C nucleotides).
        
        Returns:
            float: GC content as percentage (0-100)
        """
        if self.length == 0:
            return 0
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / self.length) * 100
    
    def get_nucleotide_count(self):
        """
        Count occurrences of each nucleotide.
        
        Returns:
            dict: Dictionary with nucleotide counts
        """
        counts = Counter(self.sequence)
        return {
            'A': counts.get('A', 0),
            'T': counts.get('T', 0),
            'G': counts.get('G', 0),
            'C': counts.get('C', 0),
            'N': counts.get('N', 0)  # Unknown nucleotides
        }
    
    def get_nucleotide_percentage(self):
        """
        Calculate percentage of each nucleotide.
        
        Returns:
            dict: Dictionary with nucleotide percentages
        """
        counts = self.get_nucleotide_count()
        percentages = {}
        for nucleotide, count in counts.items():
            percentages[nucleotide] = (count / self.length) * 100 if self.length > 0 else 0
        return percentages
    
    def reverse_complement(self):
        """
        Generate the reverse complement of the sequence.
        
        Returns:
            str: Reverse complement sequence
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(self.sequence))
    
    def get_codon_usage(self):
        """
        Calculate codon frequency (groups of 3 nucleotides).
        
        Returns:
            dict: Dictionary with codon counts
        """
        codons = Counter()
        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i+3]
            if len(codon) == 3:
                codons[codon] += 1
        return dict(codons)
    
    def print_summary(self):
        """Print a summary of sequence statistics."""
        print(f"\n{'='*50}")
        print(f"SEQUENCE ANALYSIS SUMMARY")
        print(f"{'='*50}")
        print(f"Sequence Length: {self.get_length()} bp")
        print(f"GC Content: {self.get_gc_content():.2f}%")
        print(f"\nNucleotide Composition:")
        counts = self.get_nucleotide_count()
        percentages = self.get_nucleotide_percentage()
        for nucleotide in ['A', 'T', 'G', 'C']:
            print(f"  {nucleotide}: {counts[nucleotide]} ({percentages[nucleotide]:.2f}%)")
        print(f"\nReverse Complement: {self.reverse_complement()}")
        print(f"{'='*50}\n")


def analyze_fasta_file(filename):
    """
    Analyze all sequences in a FASTA file.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        pandas.DataFrame: DataFrame with analysis results
    """
    results = []
    
    for record in SeqIO.parse(filename, "fasta"):
        analyzer = SequenceAnalyzer(str(record.seq))
        results.append({
            'ID': record.id,
            'Description': record.description,
            'Length': analyzer.get_length(),
            'GC_Content': analyzer.get_gc_content()
        })
    
    return pd.DataFrame(results)


def plot_nucleotide_composition(sequence_analyzer, title="Nucleotide Composition"):
    """
    Create a bar plot of nucleotide composition.
    
    Args:
        sequence_analyzer (SequenceAnalyzer): Analyzer object
        title (str): Title for the plot
    """
    percentages = sequence_analyzer.get_nucleotide_percentage()
    nucleotides = list(percentages.keys())
    values = list(percentages.values())
    
    plt.figure(figsize=(8, 6))
    plt.bar(nucleotides, values, color=['red', 'blue', 'green', 'yellow', 'gray'])
    plt.xlabel('Nucleotide')
    plt.ylabel('Percentage (%)')
    plt.title(title)
    plt.ylim(0, 100)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig('nucleotide_composition.png')
    print("Plot saved as 'nucleotide_composition.png'")
    plt.show()


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*50)
    print("BIOINFORMATICS ANALYSIS EXAMPLE")
    print("="*50)
    
    # Example 1: Analyze a single sequence
    print("\n1. ANALYZING A SINGLE DNA SEQUENCE:")
    dna_sequence = "ATCGATCGATCGATCGATCG"
    analyzer = SequenceAnalyzer(dna_sequence)
    analyzer.print_summary()
    
    # Example 2: Analyze multiple sequences
    print("\n2. ANALYZING MULTIPLE SEQUENCES:")
    sequences = [
        "ATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTA",
        "AAATTTGGGCCC"
    ]
    
    data = []
    for idx, seq in enumerate(sequences, 1):
        analyzer = SequenceAnalyzer(seq)
        data.append({
            'Sequence_ID': f'seq_{idx}',
            'Length': analyzer.get_length(),
            'GC_Content': analyzer.get_gc_content(),
            'A': analyzer.get_nucleotide_count()['A'],
            'T': analyzer.get_nucleotide_count()['T'],
            'G': analyzer.get_nucleotide_count()['G'],
            'C': analyzer.get_nucleotide_count()['C']
        })
    
    df = pd.DataFrame(data)
    print("\nSequence Analysis Results:")
    print(df.to_string(index=False))
    
    # Example 3: Create visualization
    print("\n3. CREATING VISUALIZATION:")
    analyzer = SequenceAnalyzer(sequences[0])
    plot_nucleotide_composition(analyzer, title="Example Sequence Composition")
    
    # Example 4: Codon usage
    print("\n4. CODON USAGE ANALYSIS:")
    codon_usage = analyzer.get_codon_usage()
    print(f"Codons found: {codon_usage}")
    
    print("\n" + "="*50)
    print("ANALYSIS COMPLETE!")
    print("="*50)
    print("\nTo use with your own data:")
    print("  1. Place FASTA files in the 'data/' directory")
    print("  2. Modify the script to load your files")
    print("  3. Run: python bio_analysis.py")
    print("="*50 + "\n")
