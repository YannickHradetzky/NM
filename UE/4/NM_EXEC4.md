# Task 1) Flops of cholesky_factorization
```python 
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
```

$$
\begin{align}
\overbrace{ \sum_{i=1}^{n} }^{ \text{first for-loop } }\overbrace{ \sum_{j=1}^{n} }^{ \text{ second for-loop} } \begin{cases}
2j  & \text{if }i=j \\
2j + 1 & \text{if }i\neq j
\end{cases} 
\end{align}
$$
We are just gonna look at the case where $i\neq j$

$$
\begin{align}
\sum_{i=1}^{n}\sum_{j=1}^{n}2j+1  & = 2\sum_{i=1}^{n}\sum_{j=1}^{n}j + \sum_{i=1}^{n}\sum_{j=1}^{n}1 \\
 & = 2n \sum_{j=1}^{n}j + n^{2} \\
 & = 2n \frac{n(n+1)}{2} + n^2 \\
 & = n^3+2n^2
\end{align}
$$
So to leading order the amount of flops scales with an upper bound of $n^3+2n^{2}$

# Task 2) Reduced QR factorization

We will compute the reduced QR factorization of the matrix 
$$
A =
\begin{bmatrix}
3 & 7 \\
0 & 12 \\
4 & 1
\end{bmatrix}
$$
using the Gram-Schmidt orthogonalization process, and then extend it to a full QR factorization.

---
1) First we normalize the first column of $A$:

$$
\begin{align}
\vec{q_{1}} = \frac{\vec{a_{1}}}{\lvert \vec{a_{1}} \rvert } = \begin{pmatrix}
\frac{3}{5} \\
0 \\
\frac{4}{5}
\end{pmatrix}
\end{align}
$$
2) Project this onto the second column of $A$:

$$
\begin{align}
\vec{v_{2}}  & = \vec{a_{2}} - (\vec{q_{1}}^{T}\cdot   \vec{a_{2}} ) \cdot  \vec{q_{1}} \\
 & =\begin{pmatrix}
7 \\
12 \\
1
\end{pmatrix} - 5 \cdot \begin{pmatrix}
\frac{3}{5} \\
0  \\
\frac{4}{5}
\end{pmatrix} = \begin{pmatrix}
4 \\
12 \\
-3
\end{pmatrix} \\
\vec{q_{2}}  & = \frac{\vec{v_{2}}}{\lvert \vec{v_{2}} \rvert } = \begin{pmatrix}
\frac{4}{13} \\
\frac{12}{13} \\
-\frac{3}{13}
\end{pmatrix}
\end{align}
$$

3) collect into Q-Matrix:
$$
\begin{align}
Q = \begin{pmatrix}
\frac{3}{5} & \frac{4}{13} \\
0 & \frac{12}{13} \\
\frac{4}{5} & -\frac{3}{13}
\end{pmatrix}
\end{align}
$$
4) Calculate the $R$ matrix:
$$
\begin{align}
r_{11} &  = \lvert \vec{a_{1}} \rvert = 5 \\
r_{12} & = \vec{q_{1}}^{T}\cdot \vec{a_{2}} = 5 \\
r_{22} & = \lvert \vec{v_{2}} \rvert = 13
\end{align}
$$
5) Build up $R$
$$
R = \begin{align}
\begin{pmatrix}
5 & 5 \\
0 & 13
\end{pmatrix}
\end{align}
$$

## Extend $Q$:
1. Compute $\vec{q_3}$, orthogonal to $\vec{q_1}$ and $\vec{q_2}$:
   $$
   \vec{q_3} = \text{normalize} \left( \vec{e_3} - (\vec{q_1}^T \cdot \vec{e_3}) \cdot \vec{q_1} - (\vec{q_2}^T \cdot \vec{e_3}) \cdot \vec{q_2} \right).
   $$

2. Collect into:
   $$
   Q_{\text{full}} =
   \begin{bmatrix}
   \frac{3}{5} & \frac{4}{13} & \cdots \\
   0 & \frac{12}{13} & \cdots \\
   \frac{4}{5} & -\frac{3}{13} & \cdots
   \end{bmatrix}.
   $$


The final result is:
$$
A = Q_{\text{full}} R.
$$

The extension was done in Python to save time!
## Python Code
```python
import numpy as np

def gram_schmidt(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
	"""
	Reduced QR factorization of a matrix A using the Gram-Schmidt orthogonalization process.
	"""
	Q = np.zeros((A.shape[0], A.shape[1]))
	R = np.zeros((A.shape[1], A.shape[1]))
	for i in range(A.shape[1]):
		Q[:, i] = A[:, i]
		for j in range(i):
			R[j, i] = np.dot(Q[:, j], A[:, i])
			Q[:, i] -= R[j, i] * Q[:, j]
		R[i, i] = np.linalg.norm(Q[:, i])
		Q[:, i] /= R[i, i]
	return Q, R

A = np.array([[3, 7], [0, 12], [4, 1]])
Q, R = gram_schmidt(A)
print("Reduced Q :\n", Q)
print("Reduced R :\n", R)

def full_qr_factorization(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
	"""
	Full QR factorization of a matrix A.
	"""
	m, n = A.shape

	# Step 1: Compute the reduced QR factorization
	Q_reduced, R = gram_schmidt(A)

	# Step 2: Extend Q to a full orthonormal basis
	Q_full = np.zeros((m, m))
	Q_full[:, :n] = Q_reduced

	# Find the orthogonal complement
	for i in range(n, m):
		# Start with a standard basis vector
		v = np.zeros(m)
		v[i] = 1

		# Subtract projections onto existing Q columns
		for j in range(i):
			v -= np.dot(Q_full[:, j], v) * Q_full[:, j]

		# Normalize the new vector
		v /= np.linalg.norm(v)
		Q_full[:, i] = v
	# Compute R_full
	R_full = Q_full.T @ A
	return Q_full, R_full

Q_full, R_full = full_qr_factorization(A)
print("Full Q:")
print(Q_full)
print("\nFull R:")
print(R_full)
```

# Task 4) backward stability of the Householder method

```python
import numpy as np
R = np.triu(np.random.rand(50,50)) 
Q = np.linalg.qr(np.random.rand(50,50))[0] 
A = Q@R

Q2, R2 = np.linalg.qr(A)

# Compute forward and backward errors
forward_error_Q = np.linalg.norm(Q - Q2) / np.linalg.norm(Q)
forward_error_R = np.linalg.norm(R - R2) / np.linalg.norm(R)
backward_error = np.linalg.norm(A - Q2@R2) / np.linalg.norm(A)

print(f"Forward error in Q: {forward_error_Q:.2e}")
print(f"Forward error in R: {forward_error_R:.2e}")
print(f"Backward error: {backward_error:.2e}")

# Compute condition number of A
cond_A = np.linalg.cond(A)
print(f"\nCondition number of A: {cond_A:.2e}")

# Interpretation
print("\nInterpretation:")
print("1. The backward error is very small, showing the computed QR")
print("   factorization is backward stable (QR ≈ A)")
print("2. The forward errors in Q and R are larger than the backward")
print("   error but still reasonable, indicating some sensitivity")
print("3. The condition number gives us the worst-case amplification")
print("   of relative errors, explaining the relationship between")
print("   forward and backward errors")
```

# Task 5) QR via Householder reflections

```python
import numpy as np

def solve_qr_householder(A, b):
    Q,R = np.linalg.qr(A)
    # if m > n, compute least squares solution
    if A.shape[0] == A.shape[1]:
        # solve Ax = b using QR factorization
        x_with_A = np.linalg.solve(A, b)
        # instead of Ax = b we use QR factorization so we solve Rx = Q^T @ b
        x_with_Q = np.linalg.solve(R, Q.T @ b)
    else:
        # solve least squares problem
        x_with_A = np.linalg.lstsq(A, b, rcond=None)[0]
        # instead of Ax = b we use QR factorization so we solve Rx = Q^T @ b
        x_with_Q = np.linalg.lstsq(R, Q.T @ b, rcond=None)[0]

    # print both solutions and the relative and absolute errors
    print(f"x_with_A: {x_with_A}")
    print(f"x_with_Q: {x_with_Q}")
    print(f"Relative error: {np.linalg.norm(x_with_A - x_with_Q) / np.linalg.norm(x_with_A)}")
    print(f"Absolute error: {np.linalg.norm(x_with_A - x_with_Q)}")
    return x_with_A, x_with_Q

A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
b = np.array([1, 2, 3])
x_with_A, x_with_Q = solve_qr_householder(A, b)
```



