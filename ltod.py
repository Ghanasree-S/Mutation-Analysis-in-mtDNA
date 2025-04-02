from collections import deque

def bfs_locus_to_disease(graph, start_node, max_depth=3):
    """
    Optimized BFS for finding diseases from a mutation (Locus → Allele → Disease).
    
    Time Complexity:
    - Best case: O(1) (if start_node is not in graph)
    - Worst case: O(B^d) ~ O(V + E) (if search explores all nodes)

    Parameters:
    - graph (dict): Graph representation of the relationships
    - start_node (str): The mutation or locus to start the search from
    - max_depth (int): The maximum depth to traverse

    Returns:
    - list: Paths leading to diseases
    """
    if start_node not in graph:
        print(f"⚠️ {start_node} not found in the graph.")
        return []

    queue = deque([[start_node]])  # Initialize queue with the starting node
    paths = []  # Stores valid disease paths
    visited = set()  # Track visited nodes

    print(f"🔎 BFS starting from: {start_node}")

    while queue:
        path = queue.popleft()  # Dequeue the first path
        node = path[-1]  # Last node in the path

        if len(path) > max_depth:  # Stop search at max depth
            continue

        if " " in node:  # Check if node is a disease (assuming diseases contain spaces)
            paths.append(path)  # Store valid path
            print(f"✅ Found Disease Path: {' → '.join(path)}")
        else:
            for neighbor in graph.get(node, []):  # Explore neighbors
                if neighbor not in visited:  # Prevent cycles
                    visited.add(neighbor)
                    new_path = path + [neighbor]  # Create new path
                    queue.append(new_path)  # Enqueue new path
                    print(f"🔄 Expanding Path: {' → '.join(new_path)}")

    return paths
