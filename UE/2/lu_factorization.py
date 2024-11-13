import numpy as np

def lu_factorization(A):
    """
    Compute the LU factorization of matrix A using Gaussian elimination with partial pivoting.
    Returns matrices P, L, and U such that PA = LU, where:
    - P is the permutation matrix
    - L is lower triangular with ones on the diagonal
    - U is upper triangular
    
    Parameters:
    A : ndarray
        Input n x n matrix
        
    Returns:
    P : ndarray
        Permutation matrix
    L : ndarray
        Lower triangular matrix
    U : ndarray
        Upper triangular matrix
    """
    n = len(A)
    
    # Initialize matrices
    U = A.copy().astype(float)
    L = np.eye(n)
    P = np.eye(n)
    
    # Gaussian elimination with partial pivoting
    for k in range(n-1):
        # Find pivot
        pivot_idx = k + np.argmax(abs(U[k:, k]))
        
        # Skip if the pivot is zero (matrix is singular)
        if abs(U[pivot_idx, k]) < 1e-10:
            raise ValueError("Matrix is singular or nearly singular")
            
        # Swap rows if necessary
        if pivot_idx != k:
            # Swap rows in U
            U[[k, pivot_idx]] = U[[pivot_idx, k]]
            # Swap rows in P
            P[[k, pivot_idx]] = P[[pivot_idx, k]]
            # Swap the part of L that's already been computed
            if k > 0:
                L[[k, pivot_idx], :k] = L[[pivot_idx, k], :k]
        
        # Compute multipliers
        for i in range(k+1, n):
            mult = U[i, k] / U[k, k]
            L[i, k] = mult
            U[i, k:] -= mult * U[k, k:]
    
    return P, L, U

def verify_factorization(A, P, L, U):
    """
    Verify that PA = LU
    Returns the maximum absolute difference between PA and LU
    """
    PA = P @ A
    LU = L @ U
    return np.max(np.abs(PA - LU))

# Example usage
if __name__ == "__main__":
    # Create a sample matrix
    A = np.array([
        [2, 1, 1],
        [4, -6, 0],
        [-2, 7, 2]
    ])
    
    # Compute LU factorization
    P, L, U = lu_factorization(A)
    
    print("Original matrix A:")
    print(A)
    print("\nPermutation matrix P:")
    print(P)
    print("\nLower triangular matrix L:")
    print(L)
    print("\nUpper triangular matrix U:")
    print(U)
    
    # Verify the factorization
    error = verify_factorization(A, P, L, U)
    print(f"\nMaximum error in PA = LU: {error}")