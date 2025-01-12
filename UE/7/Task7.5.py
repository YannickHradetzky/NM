import numpy as np

# Define the adjacency matrix for the webpage graph
# Order of pages: A, B, C, D, E
adjacency_matrix = np.array([
    [0, 1, 0, 1, 1],  # A's outgoing links
    [1, 0, 1, 0, 0],  # B's outgoing links
    [1, 1, 0, 0, 0],  # C's outgoing links
    [0, 1, 0, 0, 0],  # D's outgoing links
    [1, 0, 1, 0, 0]   # E's outgoing links
])

def calculate_pagerank(p, max_iterations=100, tolerance=1e-8):
    n = len(adjacency_matrix)
    
    # Normalize the adjacency matrix
    out_degrees = np.sum(adjacency_matrix, axis=1)
    transition_matrix = adjacency_matrix / out_degrees[:, np.newaxis]
    
    # Handle dangling nodes (if any)
    transition_matrix = np.nan_to_num(transition_matrix)
    
    # Initialize PageRank values
    pagerank = np.ones(n) / n
    
    for _ in range(max_iterations):
        prev_pagerank = pagerank.copy()
        
        # PageRank calculation with damping factor
        pagerank = (1 - p) * transition_matrix.T.dot(prev_pagerank) + p/n * np.ones(n)
        
        # Check convergence
        if np.sum(np.abs(pagerank - prev_pagerank)) < tolerance:
            break
    
    return pagerank

# Calculate PageRank for p = 0.15
p1 = 0.15
pagerank1 = calculate_pagerank(p1)
print(f"\nPageRank values with p = {p1}:")
for i, rank in enumerate(['A', 'B', 'C', 'D', 'E']):
    print(f"Page {rank}: {pagerank1[i]:.4f}")

# Calculate PageRank for p = 1e-03
p2 = 1e-03
pagerank2 = calculate_pagerank(p2)
print(f"\nPageRank values with p = {p2}:")
for i, rank in enumerate(['A', 'B', 'C', 'D', 'E']):
    print(f"Page {rank}: {pagerank2[i]:.4f}")
