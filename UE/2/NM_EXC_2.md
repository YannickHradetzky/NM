#NM 

# 1) Condition of $\cos(x)$

**Condition of $\cos(x)$**

$$
\begin{align}
\kappa_{\cos}  & = \left\lvert  x \frac{\frac{d}{dx}\cos(x)}{\cos(x)}  \right\rvert  \\
 & = \left\lvert  x \frac{-\sin(x)}{\cos(x)}  \right\rvert  \\
 & =\lvert x\tan(x) \rvert 
\end{align}
$$
so for $x \approx 0$ we use the small angle approximation $\tan(x\approx 0) \approx x$

$$
\kappa_{\cos} = \lvert x \tan(x) \rvert \approx x^{2}
$$
And for small $x$ this becomes very small

Assuming we have a **perfect implementation** I would assume that the backward evaluation is stable

Let's look at a small perturbation: 

$$
\begin{align}
\tilde{f}(x) = \text{fl}(\cos(x)) \text{ the computed value of } \cos(x)
\end{align}
$$
We want to know $\tilde{f}(x)$ for a small perturbation $\delta x$: 

$$
\begin{align}
\tilde{f}(x) & =\cos(x+\delta x) \underbrace{ \approx }_{ \text{1st Taylor} } \cos(x)-\delta x\sin(x) \\
 \delta x  & \approx \frac{\cos(x)-\tilde{f}(x)}{\sin(x)}
\end{align}
$$

So if the implementation is indeed perfect then **$\cos(x)-\tilde{f}(x) \approx 0$** and the **error we obtain from calculating backwards is minimal.**
***
# 2) regular, lower triangle band matrix 

## a) How does this matrix look like? 

$$
L = \begin{bmatrix}
\ell_{11} & 0 & 0 & \cdots & 0 & 0 & 0 \\
\ell_{21} & \ell_{22} & 0 & \cdots & 0 & 0 & 0 \\
\ell_{31} & \ell_{32} & \ell_{33} & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
\ell_{p+1,1} & \ell_{p+1,2} & \ell_{p+1,3} & \cdots & \ell_{p+1,p+1} & 0 & 0 \\
0 & \ell_{p+2,2} & \ell_{p+2,3} & \cdots & \ell_{p+2,p+1} & \ell_{p+2,p+2} & 0 \\
0 & 0 & \ell_{p+3,3} & \cdots & \ell_{p+3,p+1} & \ell_{p+3,p+2} & \ell_{p+3,p+3} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & \ell_{n,n-p} & \ell_{n,n-p+1} & \ell_{n,n}
\end{bmatrix}
$$

## b) Show that flops are linear in n for $Lx=b$

With forward substitution we solve the following in step $i$:

$$
x_{i} = \frac{1}{l_{ii}}\left( b_{i}-\sum_{j=1}^{n}l_{ij}x_{j} \right)
$$

but since all entries with $j>i$ are $0$ we can cut off the sum at $i-1$

$$
x_{i} = \frac{1}{l_{ii}}\left( b_{i}-\sum_{j=1}^{i-1}l_{ij}x_{j} \right) \text{ for }i=2,3,\dots,p+1
$$
Then we have to take into account that we don't have to consider the first $i-p$ values in the sum.

$$
x_{i}= \frac{1}{l_{ii}}\left( b_{i}-\sum ^{i-1}_{j=i-p}l_{ij}x_{j} \right) \text{for } i= p+2,p+3,\dots,n
$$

**Counting the Operations:**
Let's go step by step and first analyze the flops for the case where $i\leq p+1$

$$
x_{i} = \underbrace{ \frac{1}{l_{ii}} }_{ \text{1 division} }\left( b_{i}\underbrace{ - }_{ \text{1 subtraction} }\overbrace{ \sum_{j=1}^{i-1}l_{ij}x_{j} }^{ \text{i-1 multi and i-1-1 additions} } \right) \text{ for }i=1,2,\dots,p+1
$$

So for the **first part** we get:
$$
\sum_{i=1}^{p+1-1}1+1+(i-1)+(i-1-1) = \sum_{i}^{p}2i-1 = 2\overbrace{ \sum_{i=1}^{p}i }^{ p(p+1)/2 }-\sum_{i=1}^{p}1 = p^{2}
$$

Now lets look at the part where $\lvert i-j \rvert>p$
We can start the sum at $i-p$ because all entries before that are $0$.


$$
x_{i}= \underbrace{ \frac{1}{l_{ii}} }_{ \text{1 Division} }\left( b_{i}\underbrace{ - }_{ \text{1 Subtraction} }\overbrace{ \sum ^{i-1}_{j=i-p}l_{ij}x_{j}  }^{ \text{p-1 multi p-1-1 additions} }\right) \text{for } i= p+2,p+3,\dots,n
$$
**It appears that for computing within this region the amount of multiplication between the matrix row part and solution column part has a fixed amount of elements.**


So for the **second part** we get: 

$$
\begin{align}
\sum_{p+2}^{n}(1+1+(p-1)+(p-2)) = \sum_{p+2}^{n}2p-1  & = 2\sum_{p+2}^{n}p-(p+2-n) \\
 & =2p(p+2-n)-(p+2-n) \\
 & =2p^{2}+2p-2pn-p-2+n \\
 & =2p^{2}-2pn -2+n
\end{align}
$$
and is therefore linear in n.

**Computing total steps**
In forward substitution we need to compute the first element before starting the algorithm
$$
x_{1} = \frac{b_{1}}{l_{11}} \text{ which takes one flop}
$$
So we get a total of 
$$
p^{2}+2p^{2}-2pn-2+n +1 = 3p^{2}-2pn+n-1
$$

***
# 3) Formulas for backward substitution 

## a) Writing down general formulas

$$
U = \begin{bmatrix}
    u_{11} & u_{12} & u_{13} & \cdots & u_{1n} \\
    0 & u_{22} & u_{23} & \cdots & u_{2n} \\
    0 & 0 & u_{33} & \cdots & u_{3n} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & 0 & \cdots & u_{nn} \\
\end{bmatrix}\cdot
\begin{bmatrix}
x_{1} \\
x_{2} \\
x_{3} \\
\vdots \\
x_{n}
\end{bmatrix}=
\begin{bmatrix}
b_{1} \\
b_{2} \\
b_{3} \\
\vdots \\
b_{n}
\end{bmatrix}
$$

So we start by multiplying $\vec{x}$ with the last column of $U$ to yield the first value of the solution $x_{n}$:

$$
u_{nn}x_{n} = b_{n} \Leftrightarrow x_{n}=\frac{b_{n}}{u_{nn}}
$$

For all the other $x_{i}$:

$$
x_{i} = \frac{1}{u_{ii}} \cdot \left( b_{i}-\sum_{j=i+1}^{n}u_{ij}x_{j}   \right)
\text{ for } i=n-1,n-2,\dots,1
$$

## b) Proofing backward stability for $d=2$

So we are interested in the case: 

$$
\begin{bmatrix}
u_{11} & u_{12} \\
0 & u_{22}
\end{bmatrix} \cdot
\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} = 
\begin{bmatrix}
b_{1} \\
b_{2}
\end{bmatrix}
$$

This results in the following two equations for solving this system: 

$$
\begin{align}
x_{2} & =\frac{b_{2}}{u_{22}} \\
x_{1}  & =  \frac{b_{1} - u_{12}x_{2}}{u_{11}} = \frac{b_{1} - u_{12}\frac{b_{2}}{u_{22}}}{u_{11}}
\end{align}
$$
**These are the perfect unperturbed solutions**

Backward stability means that the computed solution is exactly the solution to a slightly perturbed version of the original problem, with the perturbation being small enough to be considered within machine precision.


To show backward stability, one needs to show that the computed solution $\tilde{\vec{x}}$ satisfies a perturbed system $(U+\Delta U)\tilde{x}=b$ . 

So we express the solutions as potentially slightly perturbed from the true $\vec{x}$ values: 

$$
\begin{align}
\tilde{x}_{1}  & = x_{1}\cdot(1+\epsilon_{1}) \\
\tilde{x}_{2} & =x_{2}\cdot(1+\epsilon_{2})
\end{align}
$$
Where $\epsilon_{i}$ represent some small floating point error in the calculation. Which can be interpreted as a slightly perturbed system. 

$$
\vec{\tilde{x}} = \vec{x}+\delta\vec{x}=\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} + 
\underbrace{ \begin{bmatrix}
\epsilon_{1}x_{1} \\
\epsilon_{2}x_{2}
\end{bmatrix} }_{ \Delta U\cdot\vec{x} }
$$

***
# 4) Write backward substitution program

```python
import numpy as np

def back_substitution_upper_triangular(U, b):
	n = U.shape[0]
	result = []
	result.append(b[-1] / U[-1, -1])
	
	for i in range(n-2, -1, -1):
		result.insert(0, (b[i] - np.dot(U[i, i+1:], result)) / U[i, i])

	print("Verification:")
	print("U * x =")
	print(np.array2string(np.dot(U, result), precision=4, separator=', '))
	print("\nOriginal b =")
	print(np.array2string(b, precision=4, separator=', '))
	print("\nDifference:")
	print(np.array2string(np.dot(U, result) - b, precision=4, separator=', '))
	return result

# Example Usage
U = np.array([[3, 2, 1], [0, 4, 2], [0, 0, 5]])
b = np.array([1, 2, 3])

back_substitution_upper_triangular(U, b)

```

# 5) Check if matrices admit LU factorization

So instead of performing Gauss Elimination we will apply Theorem 2.5 from the script and check if all leading principal submatrices are regular $\det A\neq0$

## a) 
$$A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 0 \\ 6 & 0 & 6 \end{pmatrix}$$
**1st leading principal submatrix (1x1)**
Well this is very simple and has a non zero determinant

**2nd leading principal submatrix (2x2)**

$$
\det\begin{pmatrix}
1 & 2 \\
4 & 5
\end{pmatrix} = 5-8 = -3 \neq 0
$$
**3rd leading principal submatrix (3x3)**

Here we actually have to compute the determinant of $A$. Done with python 

```python
import numpy as np

A = np.array([[1, 2, 3], [4, 5, 0], [6, 0, 6]])
det_A = np.linalg.det(A)

print(f"The determinant of A is {det_A}")

```

**So this matrix admits LU factorization**

## b) 

$$
B = \begin{pmatrix}
0 & 2 & 3 \\
4 & 5 & 0 \\
6 & 0 & 6
\end{pmatrix}
$$

This matrix **does not admit LU factorization** because the determiant of it's first leading principal submatrix is 0 

## c)

$$
C = \begin{pmatrix}
1 & 2 & 1 \\
3 & 6 & 3  \\
1 & 2 & 0 & 
\end{pmatrix}
$$

This one also **does not** admit to LU factorization because it's **2nd leading principal submatrix** has a zero determinant. 

## d) 

$$
D = \begin{pmatrix}
6 & 2 & 1 & 0 \\
2 & 5 & 1 & 1   \\
1 & 1 & 4 & 1  \\
1 & 1 & 1 & 4
\end{pmatrix}
$$


for the third and forth I am again gonna use python to solve the determiant: 

```python
import numpy as np

third = np.array([[6, 2, 1], [2, 5, 1], [1, 1, 4]])
fourth = np.array([[6, 2, 1, 0], [2, 5, 1, 1], [1, 1, 4,1], [1, 1, 1, 4]])
det_third = np.linalg.det(third)
det_fourth = np.linalg.det(fourth)

print(f"The determinant of the third principal submat = {det_third}")
print(f"The determinant of the fourth principal submat = {det_fourth}")

```
Both of them are non zero and therefore this matrix **admits** LU factorization


# 6) Solve linear Systems

## a) 
Given the system:

$$
\begin{pmatrix}
2 & -2 & 1 \\
0 & 1 & 2 \\
5 & 3 & 1
\end{pmatrix} 
\cdot 
\begin{pmatrix}
x_{1} \\
x_{2} \\
x_{3} 
\end{pmatrix} 
= 
\begin{pmatrix}
6 \\
3 \\
4
\end{pmatrix}
$$

we can represent this as the augmented matrix:

$$
\begin{pmatrix}
2 & -2 & 1 & | & 6 \\
0 & 1 & 2 & | & 3 \\
5 & 3 & 1 & | & 4
\end{pmatrix}
$$

**Make the first element of the first row equal to 1**
To simplify, divide the first row by 2:

$$
\begin{pmatrix}
1 & -1 & 0.5 & | & 3 \\
0 & 1 & 2 & | & 3 \\
5 & 3 & 1 & | & 4
\end{pmatrix}
$$

**Make all elements below the leading 1 in the first column zero**
Replace the third row with (Row 3 - 5 * Row 1):

$$
\begin{pmatrix}
1 & -1 & 0.5 & | & 3 \\
0 & 1 & 2 & | & 3 \\
0 & 8 & -1.5 & | & -11
\end{pmatrix}
$$

**Make all elements below the leading 1 in the second column zero**
Replace the third row with (Row 3 - 8 * Row 2):

$$
\begin{pmatrix}
1 & -1 & 0.5 & | & 3 \\
0 & 1 & 2 & | & 3 \\
0 & 0 & -17.5 & | & -35
\end{pmatrix}
$$

**Make the leading coefficient of the third row a 1**
Divide the third row by $-17.5$:

$$
\begin{pmatrix}
1 & -1 & 0.5 & | & 3 \\
0 & 1 & 2 & | & 3 \\
0 & 0 & 1 & | & 2
\end{pmatrix}
$$

**Back-substitute to make all elements above the leading 1s zero**
Replace Row 2 with (Row 2 - 2 * Row 3):

$$
\begin{pmatrix}
1 & -1 & 0.5 & | & 3 \\
0 & 1 & 0 & | & -1 \\
0 & 0 & 1 & | & 2
\end{pmatrix}
$$

Then replace Row 1 with (Row 1 - 0.5 * Row 3):

$$
\begin{pmatrix}
1 & -1 & 0 & | & 2 \\
0 & 1 & 0 & | & -1 \\
0 & 0 & 1 & | & 2
\end{pmatrix}
$$

Finally, replace Row 1 with (Row 1 + Row 2):

$$
\begin{pmatrix}
1 & 0 & 0 & | & 1 \\
0 & 1 & 0 & | & -1 \\
0 & 0 & 1 & | & 2
\end{pmatrix}
$$

**Solution**
Thus, the solution is:

$$
x_1 = 1, \quad x_2 = -1, \quad x_3 = 2
$$

## b) 
Given the system:

$$
\begin{pmatrix}
1 & -2 & -3 \\
5 & 6 & -1 \\
1 & -1 & -3
\end{pmatrix} 
\cdot 
\begin{pmatrix}
x_{1} \\
x_{2} \\
x_{3} 
\end{pmatrix} 
= 
\begin{pmatrix}
10 \\
2 \\
6
\end{pmatrix}
$$

we can represent this as the augmented matrix:

$$
\begin{pmatrix}
1 & -2 & -3 & | & 10 \\
5 & 6 & -1 & | & 2 \\
1 & -1 & -3 & | & 6
\end{pmatrix}
$$

**Make the first element of the first row equal to 1**
The first row already has a leading 1, so we can move to the next step.

**Make all elements below the leading 1 in the first column zero**
Replace Row 2 with (Row 2 - 5 * Row 1) and Row 3 with (Row 3 - Row 1):

$$
\begin{pmatrix}
1 & -2 & -3 & | & 10 \\
0 & 16 & 14 & | & -48 \\
0 & 1 & 0 & | & -4
\end{pmatrix}
$$

**Make the leading coefficient of the second row a 1**
Divide Row 2 by 16:

$$
\begin{pmatrix}
1 & -2 & -3 & | & 10 \\
0 & 1 & 0.875 & | & -3 \\
0 & 1 & 0 & | & -4
\end{pmatrix}
$$

**Make all elements below and above the leading 1 in the second column zero**
Replace Row 1 with (Row 1 + 2 * Row 2) and Row 3 with (Row 3 - Row 2):

$$
\begin{pmatrix}
1 & 0 & -1.25 & | & 4 \\
0 & 1 & 0.875 & | & -3 \\
0 & 0 & -0.875 & | & -1
\end{pmatrix}
$$

**Make the leading coefficient of the third row a 1**
Divide Row 3 by $-0.875$:

$$
\begin{pmatrix}
1 & 0 & -1.25 & | & 4 \\
0 & 1 & 0.875 & | & -3 \\
0 & 0 & 1 & | & \approx 1.14
\end{pmatrix}
$$

**Back-substitute to make all elements above the leading 1s zero**
Replace Row 1 with (Row 1 + 1.25 * Row 3) and Row 2 with (Row 2 - 0.875 * Row 3):

$$
\begin{pmatrix}
1 & 0 & 0 & | & 5.43 \\
0 & 1 & 0 & | & -4 \\
0 & 0 & 1 & | & 1.14
\end{pmatrix}
$$

**Solution**
Thus, the solution is:

$$
x_1 \approx 5.43, \quad x_2 = -4, \quad x_3 \approx 1.14
$$

# 7) Gaussian ellimination with partial pivoting

Not interested in that

# 8) Write program for LU factorization with partial pivoting

```python
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
```