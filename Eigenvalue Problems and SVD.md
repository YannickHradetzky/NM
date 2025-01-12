In this chapter we consider a matrix $A\in C^{n\times n}$ and numerically find its eigenvalues

# Mathematical Background
For a given matrix A we search pairs of **eigenvalues** and **eigenvectors** $x \in C^{n}$ that satisfy
$$
Ax=\lambda x
$$
- $x \in C^{n} \neq 0$
- $\lambda \in C$

Analytically we find the **eigenvalues** by the roots $\lambda_{i}$ of the characteristic polynomial:
$$
\chi(z)=\det(A-zI)
$$
which is a non-linear problem. 
- **Multiplicity** of an eigenvalue $\lambda$ is its multiplicity as a root of $\chi(z)$

For each eigenvalue $\lambda$ we can calculate the corresponding eigenvectors by solving the linear equation
$$
(A-\lambda I)v_{j}=0
$$

The linear space **spanned** by the eigenvectors $v_{j}$ is called the **eigenspace** and is denoted with $E_{\lambda}$ . The number of linear independent eigenvectors (that is the dimension of $E_{\lambda}$) is called the geometric multiplicity of the eigenvalue $\lambda$.

- **positive definite** - All eigenvalues are **real** and **positive**
- **hermitian** ($A=A^{*}$ ) - all eigenvalues are **real**
- **non singular** - If $x$ is an eigenvector then $\frac{1}{\lambda}$ is an eigenvalue of $A^{-1}$

---
# Estimation of eigenvalues

For a quick estimation one can use [[Gershgorin circle theorem]] to get a realm of where the eigenvalues might be.  To further confine the range in which the eigenvalues lie once can use the [[Rayleigh quotient]]

## Lemma 5.2
Let $A$ be hermitian and $r(x)$ the **Rayleigh quotient**. Then the eigenvectors $v_{1},v_{2},\dots v_{n}$ of $A$ are the stationary points of $r(x)$ and the corresponding eigenvalues $\lambda_{1},\lambda_{2},\dots,\lambda_{n}$ of $A$ are the stationary values. Further there holds
$$
r(x)-r(v_{i})=O(\lvert \lvert v_{i}-x \rvert  \rvert ^{2})
$$
as $x\to v_{i}$

## The range of $A$
$$
W(A):=\{ r(x) : x \in C^{n} \}
$$
- is connected
- If $A$ is hermitian, then the range of $A$ is the real interval $[\lambda_{min},\lambda_{max}]$
- If $A$ is **skew-symmetric** ($A^{*}=-A$), then the range of $A$ is an imaginary interval.

Every matrix $A \in C^{n\times n}$ can be split into an hermitian part and a skew-symmetric part
$$
A = \underbrace{ \frac{A+A^{*}}{2} }_{ \text{hermitian} }+\underbrace{ \frac{A-A^{*}}{2} }_{ \text{skew-symmetric} }
$$
Using this splitting and the above properties we directly derive: 

---
## Theorem of Bendixon
$$
\sigma(A) \subseteq W\left( \frac{A+A^{*}}{2} \right) \oplus W\left( \frac{A-A^{*}}{2} \right)
$$
where $\sigma(A)$ is the spectrum of $A \in C^{n\times n}$.

---

# Numerical treatment of eigenvalue problems

## Power iteration 

How to numerically calculate eigenvalues and eigenvectors by a method called **power iteration**. 
We only consider symmetric matrices and sort the **absolute value** of  their eigenvalues in **descending order**.

### Algorithm 1 (Power Iteration)
Approximate the eigenvector of a matrix corresponding to the largest eigenvalue
Approximation of eigenvalue is possible with [[Rayleigh quotient]]

```python
A = np.rand(n,n)
sym_mat_A = (A + A.T)/2
A = sym_mat_A

z = np.rand(z) # init z vector

while converged != True:
	new_z = np.dot(A,z)
	new_z /= np.linalg.norm(new_z)
	z = new_z
```
---
## Inverse and Rayleigh quotient iteration 
### Algorithm 2 (Inverse iteration)
Is an improvement of the **power iteration** algorithm by utilizing that 
$$
\begin{align}
A \vec{v} & = \lambda \vec{v} \\
(A-\mu I)\vec{v} & =(\lambda-\mu)\vec{v} \\
(\lambda-\mu)^{-1}\vec{v} & =(A-\mu I)^{-1}\vec{v}
\end{align}
$$
If we chose a value $\mu$ that is close to an eigenvalue $\lambda$. Then $(\lambda-\mu)^{-1}$ will be large, and therefore the eigenvector of $(A-\mu I)^{-1}$ can be computed very fast by power iteration. 
It can be further improved by successively improving the eigenvalue estimate $\mu$ in each step. 
```python
matrix_A = np.rand(n,n)
A = (matrix_A + matrix_A.T)/2

z_0 = np.rand(n)
z_0 /= np.linalg.norm(z_0)

while converged != True:
	# Solve (A-mu*I)w = z_k for w
	z_k = w/np.linalg.norm(w)
```
Further this can be improved by calculating the [[Rayleigh quotient]] in every step (of vector $w$) and use it as a new approximation. 

--- 
## QR algorithm 
Is a stable and simple procedure for calculating all eigenvalues and eigenvectors. 

### Algorithm
```python
Input: Matrix A
Output: Matrix A with eigenvalues in the diagonal and Q as an orthogonal matrix with eigenvectors as columns

# Initialization
k = 0
A_k = A
Q_total = I  # Identity matrix, to accumulate orthogonal transformations

# Iterate until convergence criterion is satisfied
while not ConvergenceCriterion(A_k) do
    # Perform QR factorization
    Q_k, R_k = QRFactorization(A_k)
    # Update A_k for the next iteration
    A_k = R_k * Q_k  # Recombine in reverse order
    # Accumulate orthogonal transformations
    Q_total = Q_total * Q_k
    # Increment iteration counter
    k = k + 1
end

# Final results
A = A_k  # A contains eigenvalues on the diagonal
Q = Q_total  # Q contains eigenvectors as columns
```

- For a **hermitian** matrix the algorithm converges to a diagonal matrix of eigenvalues. 
- The matrix $Q=Q_{0}\dots Q_{k}$ is **orthogonal** and the columns are **approximations** of the eigenvectors.
- For the **non-hermitian** case, the algorithm converges to the (real) [[Schur-form]] of the matrix. 

**Enhancements**
- Transform $A$ to its **tridiagonalform** using the **Householder transform**.
- Instead of $A_{k}$, the shifted matrix $A_{k}-\mu_{k}I$ is factored, where $\mu_{k}$ is an eigenvalue estimate.
	- Similar idea as Rayleigh quotient iteration. 
- If possible $A_{k}$ is split into sub-matrices

---
## Deflation and dimension reduction

Let's again assume a **symmetric matrix $A$.**  Once could compute the eigenvector $\vec{v_{1}}$ (normalized) and it's corrisponding (largest absolute eigenvalue) $\lambda_{1}$ by power iteration and apply deflation like:
$$
\tilde{A}=A-\lambda_{1}\vec{v_{1}}\vec{v_{1}}^{*}
$$
- $\tilde{A}$ has the same eigenvalues and eigenvectors but the largest is shifted to zero.
- Now all other eigenvectors can be computed using $\tilde{A}$. 


# Generalized eigenvalue problem 

## Definition:
The problem of finding a vector $\vec{v}\in C^{n}/ \{ 0 \} \}$ so that $A \vec{v}=\lambda B \vec{v}$  for $A,B \in C^{n\times n}$ is called a **generalized eigenvalue problem** and the **generalized eigenvalue** $\lambda \in C$ obeys:
$$
\det(A-\lambda_{i}B)=0
$$
If $B$ is non-singular, the solution can be obtained by solving the equivalent standard eigenvalue problem:

$$
B^{-1}A \vec{v}=\lambda  \vec{v}
$$

If $A,B$ are real and symmetric and $B$ is also positive definite, we can define a generalization to the [[Rayleigh quotient]]:
$$
r(x):=\frac{x^{*}Ax}{x^{*}Bx}
$$
from $\nabla r(x)=0$ follows that $Ax=r(x)B$.

Thus, the extreme points of the generalized **Rayleigh quotient** are attained by the eigenvectors of the generalized eigenvalue problems. 

With 
- $\Lambda$ a diagonal matrix containing the eigenvalues
- $V$ a matrix with the eigenvalues as columns
one can rewrite the generalized eigenvalue problem as: 
$$
AV=BV\Lambda
$$

# Singular Value Decomposition (SVD)

Singular Value Decomposition (SVD) is a powerful matrix factorization technique used in linear algebra, data compression, and machine learning. It decomposes a given matrix into three simpler matrices, revealing essential properties of the matrix.

## Definition
For a given matrix $A$ of dimensions $m \times n$, the Singular Value Decomposition is:

$$A = U \Sigma V^T$$

Where:
- $U$ is an $m \times m$ **orthogonal matrix** (columns are orthonormal vectors).
- $\Sigma$ is an $m \times n$ **diagonal matrix** containing the singular values of $A$ (non-negative real numbers) in descending order.
- $V^T$ (or $V^\dagger$ for complex matrices) is the transpose (or conjugate transpose) of an $n \times n$ **orthogonal matrix**.

## Components Breakdown
1. **Singular Values ($\Sigma$):** The singular values are the square roots of the eigenvalues of $A^T A$ or $A A^T$. They represent the magnitude of the directions of transformation induced by $A$.
2. **Left Singular Vectors ($U$):** Columns of $U$ are eigenvectors of $A A^T$, representing the directions in the input space that $A$ stretches or compresses.
3. **Right Singular Vectors ($V$):** Columns of $V$ are eigenvectors of $A^T A$, corresponding to the directions in the output space.

## Steps to Compute SVD
1. Compute $A^T A$ and find its eigenvalues and eigenvectors:
   - Eigenvalues of $A^T A$ are the squared singular values ($\sigma^2$) of $A$.
   - Eigenvectors of $A^T A$ form the columns of $V$.
2. Compute $A A^T$ and find its eigenvalues and eigenvectors:
   - Eigenvalues match those of $A^T A$.
   - Eigenvectors of $A A^T$ form the columns of $U$.
3. Arrange the singular values (square roots of eigenvalues) in descending order to form $\Sigma$.

## Properties of SVD
1. **Low-Rank Approximation:** SVD can approximate a matrix by keeping only the largest $k$ singular values and their corresponding vectors. This is foundational for techniques like Principal Component Analysis (PCA) and Latent Semantic Analysis (LSA).
2. **Compression:** SVD enables compression by storing only the largest singular values and their corresponding singular vectors.
3. **Stability:** SVD is numerically stable and often used to compute pseudoinverses and solve least squares problems.
4. **Energy Capture:** The sum of squares of singular values represents the total energy (variance) in the data. Retaining top $k$ singular values captures most of the matrix's information.

## Geometric Interpretation
- $U$: Defines a new set of orthonormal basis vectors in the input space.
- $\Sigma$: Scales or stretches the space along those basis vectors.
- $V^T$: Rotates the space in the output domain.

## Example
Let:
$$A = \begin{bmatrix} 3 & 1 \\ 1 & 3 \end{bmatrix}$$
1. Compute $A^T A$:
$$A^T A = \begin{bmatrix} 10 & 4 \\ 4 & 10 \end{bmatrix}$$
2. Eigenvalues of $A^T A$: $\lambda_1 = 14$, $\lambda_2 = 6$. Singular values: $\sigma_1 = \sqrt{14}$, $\sigma_2 = \sqrt{6}$.
3. Compute eigenvectors for $A^T A$ and $A A^T$ to form $U$ and $V$.
4. Assemble $\Sigma$ with $\sqrt{14}$ and $\sqrt{6}$ along the diagonal:
$$\Sigma = \begin{bmatrix} \sqrt{14} & 0 \\ 0 & \sqrt{6} \end{bmatrix}$$

## Applications
1. **Dimensionality Reduction:** Reduce the size of data while retaining its essential features.
2. **Data Compression:** Store only the significant components of a matrix.
3. **Noise Reduction:** Filter out low-energy components (small singular values).
4. **Image Compression:** Approximate image matrices using top singular values.
5. **Recommender Systems:** SVD is used to identify latent features in collaborative filtering models.

## Theorems and Remarks

- If $A$ is **hermitian** the singular values of $A$ are the absolute values of the eigenvalues of $A$
- If **$A$ is hermitian and positive (semi) definite** then the SVD coincides with its spectral decomposition $A=Q\Lambda Q$ 
- Stable SVD algorithms amount to an eigenvalue decomposition of $M$ where $M$:
$$
M = \begin{pmatrix}
0 & A^{*} \\
A & 0
\end{pmatrix}
$$
