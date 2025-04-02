from flask import Flask, render_template, request
from needleman_wunsch import global_alignment
from smith_waterman import local_alignment
from gene_labeling import label_gene_regions
from graph_module import Graph, bfs_locus_to_disease
import json

app = Flask(__name__)

# Load reference wild-type sequence
with open("originalmtdna/original_mtdna.fasta", "r") as f:
    wild_seq = "".join(line.strip() for line in f if not line.startswith(">"))

print(f"✅ Wild-type sequence loaded. First 50 bases: {wild_seq[:50]}...\n")

# Load reference gene locus positions
with open("mt_locus.json") as f:
    gene_locus = json.load(f)

# Initialize the Graph for disease prediction
graph = Graph()
graph.load_from_csv("mito.csv")  # Ensure mito.csv is present in the same directory

@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        file = request.files.get("usersequence")
        if not file:
            return "Error: No file uploaded", 400

        mutated_seq = file.read().decode("utf-8").strip()
        mutated_seq = "".join(line.strip() for line in mutated_seq.split("\n") if not line.startswith(">"))

        print(f"✅ Mutated sequence uploaded. First 50 bases: {mutated_seq[:50]}...\n")
        print("🔄 Performing Needleman-Wunsch Global Alignment...\n")

        # Perform global alignment
        wild_aligned, mutated_aligned, global_score = global_alignment(wild_seq, mutated_seq)
        print(f"✅ Global Alignment Complete. Score: {global_score}\n")

        # Identify gene regions
        labeled_sequences = label_gene_regions(wild_aligned, mutated_aligned)

        # Perform local alignment (Smith-Waterman) for all loci
        local_results = {}
        mutation_impact = {}

        print("🔄 Performing Smith-Waterman Local Alignment for all loci...\n")
        for locus, sequences in labeled_sequences.items():
            wild_locus_seq = sequences["wild_seq"]
            mutated_locus_seq = sequences["mutated_seq"]

            local_aligned_seq1, local_aligned_seq2, local_score = local_alignment(wild_locus_seq, mutated_locus_seq)

            # Detect mutations and format them correctly
            mutations = []
            for i, (a, b) in enumerate(zip(local_aligned_seq1, local_aligned_seq2)):
                if a != b:
                    position = gene_locus[locus][0] + i  # Adjust position to match genome
                    if a == "-":  # Insertion
                        mutation = f"m.{position}_{position+1}ins{b}"
                    elif b == "-":  # Deletion
                        mutation = f"m.{position}{a}_del"
                    else:  # Substitution
                        mutation = f"m.{position}{a}>{b}"
                    mutations.append(mutation)

            num_mutations = len(mutations)
            mutation_impact[locus] = {"score": local_score, "mutations": num_mutations, "alleles": mutations}

            local_results[locus] = {
                "score": local_score,
                "aligned_seq1": local_aligned_seq1,
                "aligned_seq2": local_aligned_seq2,
                "mutations": mutations
            }

            print(f"✅ Local Alignment Complete for {locus}. Score: {local_score}, Mutations: {num_mutations}")
            print(f"🧬 Detected Mutations: {mutations}\n")

        # Select the most impactful locus (lowest score + highest mutations)
        ranked_loci = sorted(mutation_impact.keys(), key=lambda x: (mutation_impact[x]['score'], -mutation_impact[x]['mutations']))
        top_locus = ranked_loci[0] if ranked_loci else None
        top_allele = None
        predicted_diseases = set()

        if top_locus:
            # Select the most impactful allele within this locus
            alleles = mutation_impact[top_locus]["alleles"]
            if alleles:
                top_allele = alleles[0]  # Select the first most impactful allele

                print(f"🚨 Most impactful locus identified: {top_locus}")
                print(f"🔬 Most impactful allele: {top_allele}")

                # Check if the mutation exists in the graph
                if top_allele in graph.graph:
                    print(f"✅ {top_allele} found in graph. Running BFS for disease prediction...\n")
                    
                    # Predict diseases using BFS traversal
                    disease_paths = bfs_locus_to_disease(graph.graph, top_allele)
                    print(f"🦠 BFS Disease Paths: {disease_paths}")

                    # Extract disease names from paths
                    for path in disease_paths:
                        if len(path) > 1:  # Ensure it's a valid path
                            predicted_diseases.add(path[-1])  # Last element should be the disease
                else:
                    print(f"❌ {top_allele} NOT found in graph. No disease predictions possible.\n")

        # Debugging: Ensure diseases are predicted
        print(f"⚠️ Predicted Diseases: {predicted_diseases}")

        return render_template("result.html",
                               global_score=global_score,
                               aligned_seq1=wild_aligned,
                               aligned_seq2=mutated_aligned,
                               local_results=local_results,
                               mutation_impact=mutation_impact,
                               top_locus=top_locus,
                               top_allele=top_allele,
                               predicted_diseases=list(predicted_diseases))

    return render_template("index.html")

if __name__ == "__main__":
    print("🚀 Flask server is running at http://127.0.0.1:5000/\n")
    app.run(debug=True)
