# Codon table for translation
codon_table = {
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
    'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine', 'AUG': 'Methionine',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
}

# Function to transcribe DNA to RNA
def transcribe_dna_to_rna(dna_seq):
    return dna_seq.replace('T', 'U')

# Function to translate RNA to protein
def translate_rna_to_protein(rna_seq):
    return [codon_table[rna_seq[i:i+3]] for i in range(0, len(rna_seq), 3) if codon_table[rna_seq[i:i+3]] != 'STOP']

# Function to calculate GC content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

# Function to calculate melting temperature (Tm) for DNA sequences
def calculate_tm(sequence):
    at_count = sequence.count('A') + sequence.count('T')
    gc_count = sequence.count('G') + sequence.count('C')
    return (2 * at_count) + (4 * gc_count)

# Function to design a primer for a given DNA sequence
def design_primer(sequence, length=20, gc_min=40, gc_max=60):
    """
    Designs a primer based on the given sequence.
    Parameters:
    - sequence: DNA sequence
    - length: Desired length of the primer (default 20 bases)
    - gc_min: Minimum GC content (default 40%)
    - gc_max: Maximum GC content (default 60%)
    
    Returns:
    - primer sequence or None if no suitable primer is found
    """
    for i in range(len(sequence) - length + 1):
        primer = sequence[i:i+length]
        gc_content = calculate_gc_content(primer)
        if gc_min <= gc_content <= gc_max:
            return primer
    return None  # No suitable primer found

# Example Usage
if __name__ == "__main__":
    dna_sequence = "ATGCGATAGCTAG"
    rna_sequence = transcribe_dna_to_rna(dna_sequence)
    protein_sequence = translate_rna_to_protein(rna_sequence)
    gc_content = calculate_gc_content(dna_sequence)
    melting_temp = calculate_tm(dna_sequence)

    print(f"Original DNA Sequence: {dna_sequence}")
    print(f"Transcribed RNA Sequence: {rna_sequence}")
    print(f"Protein Sequence: {', '.join(protein_sequence)}")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"Melting Temperature (Tm): {melting_temp}°C")

    # Primer Design Example
    primer = design_primer(dna_sequence)
    if primer:
        print(f"Designed Primer: {primer}")
    else:
        print("No suitable primer found.")
