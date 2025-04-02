import csv
from collections import defaultdict, deque

class Graph:
    def __init__(self):
        """Initialize an empty adjacency list graph."""
        self.graph = defaultdict(set)

    def add_edge(self, node1, node2):
        """Adds a bidirectional edge between node1 and node2."""
        self.graph[node1].add(node2)
        self.graph[node2].add(node1)

    def load_from_csv(self, file_path):
        """Builds the graph from a CSV file where:
        - Loci are connected to alleles.
        - Alleles are connected to diseases.
        """
        with open(file_path, "r") as file:
            reader = csv.reader(file)
            next(reader)  # Skip header row

            for row in reader:
                if len(row) < 3:
                    continue  # Skip incomplete rows

                locus, allele, disease = row[:3]  # Extract relevant columns
                self.add_edge(locus.strip(), allele.strip())  # Locus → Allele
                self.add_edge(allele.strip(), disease.strip())  # Allele → Disease

    def get_neighbors(self, node):
        """Returns neighbors of a given node in the graph."""
        neighbors = self.graph.get(node, set())
        print(f"📊 Graph Neighbors of {node}: {neighbors}")  # Debugging output
        return neighbors

def bfs_locus_to_disease(graph, start_node, max_depth=3):
    """
    Optimized BFS for Locus → Allele → Disease with depth limit.
    """
    if start_node not in graph:
        print(f"❌ {start_node} not found in the graph.")
        return []

    queue = deque([[start_node]])  # Initialize queue
    paths = []  # Store valid disease paths
    visited = set()  # Track visited nodes

    print(f"🔎 BFS starting from: {start_node}")

    while queue:
        path = queue.popleft()  # Dequeue the first path
        node = path[-1]  # Last node in the path

        if len(path) > max_depth:
            continue

        if " " in node:  # Check if it's a disease (assuming diseases contain spaces)
            paths.append(path)
            print(f"✅ Found Disease Path: {' → '.join(path)}")
        else:
            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = path + [neighbor]
                    queue.append(new_path)
                    print(f"🔄 Expanding Path: {' → '.join(new_path)}")

    return paths
