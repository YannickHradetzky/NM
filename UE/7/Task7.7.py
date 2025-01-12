import numpy as np

def qr_algorithm(A, max_iter=100, tol=1e-10):
    """
    Basic QR algorithm for computing eigenvalues of a matrix
    
    Parameters:
    A (ndarray): Input matrix (should be square)
    max_iter (int): Maximum number of iterations
    tol (float): Tolerance for convergence
    
    Returns:
    ndarray: Approximate eigenvalues
    """
    n = A.shape[0]
    A_k = A.copy()
    
    for k in range(max_iter):
        # Store previous diagonal for convergence check
        prev_diag = np.diag(A_k).copy()
        
        # QR decomposition
        Q, R = np.linalg.qr(A_k)
        
        # Update A_k = RQ
        A_k = R @ Q
        
        # Check convergence
        if np.allclose(prev_diag, np.diag(A_k), rtol=tol):
            break
    
    return np.diag(A_k)

# Generate a random matrix
n = 4
random_matrix = np.random.rand(n, n)

# Make it symmetric: A = (A + A^T)/2
A = (random_matrix + random_matrix.T) / 2

print("Original symmetric matrix:")
print(A)
print("\n")

# Compute eigenvalues using our QR algorithm
qr_eigenvalues = np.sort(qr_algorithm(A))
print("Eigenvalues from QR algorithm:")
print(qr_eigenvalues)
print("\n")

# Compare with numpy's eigenvalue solver
numpy_eigenvalues = np.sort(np.linalg.eigvals(A))
print("Eigenvalues from numpy.linalg.eigvals:")
print(numpy_eigenvalues)
print("\n")

# Compare the results
print("Absolute difference between methods:")
print(np.abs(qr_eigenvalues - numpy_eigenvalues))
