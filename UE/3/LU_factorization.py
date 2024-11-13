import numpy as np

def LU_factorization(A: np.ndarray) -> np.ndarray:
    """
    LU factorization of a square matrix A. 
    Stores both L and U in the same matrix, where:
    - Diagonal and upper triangle contains U
    - Lower triangle contains the multipliers of L (excluding diagonal)
    """
    n = A.shape[0]
    # Make a copy to avoid modifying input
    LU = A.copy()

    for k in range(0,n-1):
        if A[k, k] == 0:
            raise ValueError("Matrix is singular, which would lead to a division by zero.")
        
        # instead of inner looping over j we use the entire values k+1:n
        LU[k+1:,k] = LU[k+1:,k] / LU[k,k]
        LU[k+1:,k+1:] -= np.outer(LU[k+1:,k], LU[k,k+1:])  

    print(LU)   

    # extract L by taking all elements below the diagonal and adding the identity matrix
    L = np.tril(LU, k=-1) + np.eye(n)
    # extract U by taking all elements above the diagonal including the diagonal
    U = np.triu(LU)


    return L, U

if __name__ == "__main__":
    A = np.array([
        [2.0, 1.0, 1.0],
        [4.0, 3.0, 3.0],
        [8.0, 7.0, 9.0]
    ])
    L, U = LU_factorization(A)

    print(L@U==A)