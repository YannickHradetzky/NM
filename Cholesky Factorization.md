Self adjoint, positive definite matrices (SPD) always have LU factorizations. But they allow for a more simplified factorization of the form $A=R^{*}R$, where $R$ is upper triangular. **So only one matrix needs to be constructed**. 

**self-adjoint** : $A^{*}=A$
**positive definite** : All diagonal entries of SPD matrices are positive real numbers ($a_{jj}=\vec{e_{j}}^{*}A \vec{e_{j}}>0$)

[[Skript Numerics 1.pdf]] Lemma 2.3) states among other things: 
- An $m\times m$ submatrix of A is called principal submatrix of $A$, if it is obtained from $A$ by deleting any $n-m$ columns and the same $n-m$ rows. Every principal submatrix of $A$ can be written as a matrix product $J^{*}AJ$, where $J$ is obtained from the n-dimensional identity matrix by deleting some columns. Also $J^{*}AJ$ is SPD if $A$ is. 

```python
import numpy as np

def cholesky_factorization(A):
    n = A.shape[0]
    
    # Initialize the lower triangular matrix L with zeros
    L = np.zeros_like(A)

    # Print the original matrix
    print("Original matrix:")
    print(A)

    for i in range(n):
        for j in range(i + 1):
            # Compute the diagonal elements
            if i == j:
                L[i, j] = np.sqrt(A[i, i] - np.sum(L[i, :j]**2))
                print(f"Calculating L[{i},{j}] (diagonal): L[{i},{j}] = sqrt({A[i,i]} - sum(L[{i},:{j}]^2) = {L[i,j]}")
            else:
                # Compute the off-diagonal elements
                L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j])) / L[j, j]
                print(f"Calculating L[{i},{j}] (off-diagonal): L[{i},{j}] = ({A[i,j]} - sum({L[i,:j]} * {L[j,:j]}) / {L[j,j]} = {L[i,j]}")
    
    return L

# Example symmetric positive definite matrix
A = np.array([
    [4, 2, 2],
    [2, 3, 1],
    [2, 1, 3]
], dtype=float)

# Perform Cholesky factorization
L = cholesky_factorization(A)

print("\nCholesky factor L:")
print(L)
print("\nL * L.T (reconstructed matrix):")
print(np.dot(L, L.T))
```

