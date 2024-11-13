import numpy as np

def cholesky_factorization(A: np.ndarray) -> np.ndarray:
    """
    Cholesky factorization of a positive definite matrix A.
    """
    # Test square
    if A.shape[0] != A.shape[1]:
        raise ValueError("Matrix A is not square")
        
    # Test self adjoint
    if not np.allclose(A, A.T):
        raise ValueError("Matrix A is not self adjoint")
    
    # check for positive definiteness via eigenvalues
    if np.any(np.linalg.eigvals(A) <= 0):
        raise ValueError("Matrix A is not positive definite")
    
    R = np.zeros((A.shape[0], A.shape[1]))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i == j: 
                R[i, j] = np.sqrt(A[i, j] - np.sum(R[i, :j]**2))
            elif i > j:
                R[i, j] = (A[i, j] - np.sum(R[i, :j] * R[j, :j])) / R[j, j]
    return R


if __name__ == "__main__":
    # create a 5x5 SPD matrix
    A = np.array([[2,1,0,0],
                 [1,2,1,0],
                 [0,1,2,1],
                 [0,0,1,2]])
    L = cholesky_factorization(A)
    print(L)
    print(np.allclose(L @ L.T, A))