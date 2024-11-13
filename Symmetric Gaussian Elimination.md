#NM
Symmetric Gaussian Elimination Example

We’ll work through the process on a symmetric $3 \times 3$ matrix:


$$

A = \begin{bmatrix}
4 & 2 & 2 \\
2 & 3 & 1 \\
2 & 1 & 3
\end{bmatrix}
$$

1. **First Pivot**: The first pivot element is $a_{11} = 4$. We eliminate $a_{21}$ and $a_{31}$ (below the pivot) to make them zero, while ensuring that the matrix remains symmetric.

• Compute the multipliers:

$$
m_{21} = \frac{a_{21}}{a_{11}} = \frac{2}{4} = 0.5, \quad m_{31} = \frac{a_{31}}{a_{11}} = 0.5
$$
• Update the matrix:

For the second row and column:
$$
a_{22} = a_{22} - m_{21} \cdot a_{12} = 3 - 0.5 \cdot 2 = 2
$$
$$
a_{23} = a_{23} - m_{21} \cdot a_{13} = 1 - 0.5 \cdot 2 = 0
$$

For the third row and column:
$$
a_{33} = a_{33} - m_{31} \cdot a_{13} = 3 - 0.5 \cdot 2 = 2
$$
$$
a_{32} = a_{32} - m_{31} \cdot a_{12} = 1 - 0.5 \cdot 2 = 0
$$
The updated matrix is now:
$$
A = \begin{bmatrix}
4 & 2 & 2 \\
0 & 2 & 0 \\
0 & 0 & 2
\end{bmatrix}
$$

2. **Second Pivot**: Move to the next pivot element, $a_{22} = 2$.

Since the entries below this pivot are already zero (thanks to symmetry), no further elimination steps are required here.

3. **Third Pivot**: The third pivot element, $a_{33} = 2$, also does not require further elimination, as the entries below are zero.

**Resulting Triangular Form**

After applying symmetric Gaussian elimination, we’ve obtained an upper triangular matrix:

  

$$
A = \begin{bmatrix}
4 & 2 & 2 \\
0 & 2 & 0 \\
0 & 0 & 2
\end{bmatrix}

$$
This matrix is still symmetric and has been reduced to an upper triangular form through symmetric elimination steps.

```python
import numpy as np

def symmetric_gaussian_elimination(A):
    n = A.shape[0]
    
    # Print the original matrix
    print("Original matrix:")
    print(A)
    
    for i in range(n):
        # Pivot element
        pivot = A[i, i]
        print(f"\nPivoting on row {i} (pivot = {pivot}):")
        
        for j in range(i + 1, n):
            # Calculate multiplier
            if A[j, i] != 0:
                m = A[j, i] / pivot
                print(f"Eliminating element A[{j},{i}] using multiplier m = {m}")
                
                # Update row j
                for k in range(i, n):
                    A[j, k] -= m * A[i, k]
                    
                # Update the symmetric entry
                for k in range(i + 1, n):
                    A[k, j] -= m * A[k, i]

            print(f"Matrix after eliminating A[{j},{i}]:")
            print(A)
    
    return A

# Example symmetric matrix
A = np.array([
    [4, 2, 2],
    [2, 3, 1],
    [2, 1, 3]
], dtype=float)

# Perform symmetric Gaussian elimination
result = symmetric_gaussian_elimination(A)

print("\nFinal upper triangular matrix:")
print(result)
```
