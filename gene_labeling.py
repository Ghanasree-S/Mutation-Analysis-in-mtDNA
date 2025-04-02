import json

# Load mitochondrial gene reference table from JSON
with open("mt_locus.json") as f:
    gene_locus = json.load(f)

def label_gene_regions(wild_aligned, mutated_aligned):
    """
    Extract and label gene regions from the aligned sequences.
    
    Parameters:
    - wild_aligned: Aligned wild-type sequence
    - mutated_aligned: Aligned mutated sequence
    
    Returns:
    - labeled_sequences: Dictionary with gene regions, sequences, and gap counts
    """
    labeled_sequences = {}

    for gene, (start, end) in gene_locus.items():
        # Extract the portion of the aligned sequences corresponding to this gene
        wild_gene_seq = wild_aligned[start-1:end]  # Extract from aligned wild sequence
        mutated_gene_seq = mutated_aligned[start-1:end]  # Extract from aligned mutated sequence

        # Count gaps in the mutated sequence
        num_gaps = mutated_gene_seq.count('-')

        # Store the labeled sequence info
        labeled_sequences[gene] = {
            "wild_seq": wild_gene_seq,
            "mutated_seq": mutated_gene_seq,
            "gaps": num_gaps
        }

    return labeled_sequences
