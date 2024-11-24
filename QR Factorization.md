#NM 

Here we do not restrict to square matrixes. 
We consider a $A \in K^{m\times n}$ , with $m\geq n$. For such a matrix there are two different types of QR factorizations, a reduced and a full one. 

**The reduces QR factorization**
$$
A = \hat{Q}\hat{R}
$$
Here 
- $\hat{Q}$ is $m\times n$ with orthonormal columns
- $\hat{R}$ is $n\times n$ an upper triangle matrix

**The full QR factorization**
$$
A = QR
$$
Here 
- $Q \in K^{m\times m}$ is orthonormal/unitary 
- $R\in K^{m\times n}$ is generalized upper triangle
	- not necessarily square but still $r_{ij}=0 ;\forall i>j$

The relation between these two can be visualized as follows: 
# Gram Schmidt Process

QR factorization can be understood in terms of a well-known algorithm from linear algebra : **The Gram Schmidt Process**

Write down QR-factorization using column vectors:

$$
\begin{bmatrix} \\

\vec{a_{1}} & \dots & \vec{a_{n}} \\
 \\

\end{bmatrix} = 
\begin{bmatrix}
 \\
\vec{q_{1}}  & \dots & \vec{q_{n}} \\
 \\

\end{bmatrix}
\begin{bmatrix}
r_{11} & \dots & r_{1n} \\
 & \ddots &  \vdots \\
 &  & r_{nn}
\end{bmatrix}
$$

This yields the following equations :
$$
\begin{align}
\vec{a_{1}}  & = r_{11}q_{1} \\
\vec{a_{2}}  & = r_{12}q_{1}+ r_{22}q_{2} \\
 & \vdots \\
\vec{a_{n}}  & = r_{1n}q_{1} + r_{2n}q_{2}+ \dots + r_{nn}q_{n}
\end{align}
$$
So here every **column of $A$** is expressed as a **linear combination** of certain **columns of $\hat{Q}$**. Therefore the columns of $A \in \text{span}(\vec{q_{1}},\dots,\vec{q_{k}})$. Where $k$ can go till $n$.

If $A$ has full rank then also $\hat{R}$ has full rank and therefore all diagonal entries $r_{kk}$ are nonzero. 
In particular the $k\times k$ leading principal submatrix of $\hat{R}$ is invertible. So we can express $q_{k}$ as linear combinations of the first $k$ columns of $A$.  $q_{k}\in \text{span}(\vec{a_{1}},\dots,\vec{a_{k}})$ which means that the columns of $A$ and $\hat{Q}$ span the same spaces

$$
\text{span}(\vec{a_{1}},\dots,\vec{a_{k}}) = \text{span}(\vec{q_{1}},\dots,\vec{q_{k}})
$$
So we can **reformulate the problem**:
Given $n$ linearly independent vectors $\{ \vec{a_{i} } \}_{i=1}^{n}$ find $n$ orthonormal vectors $\{ \vec{q_{i}} \}_{i=1}^{n}$ spanning the same space. (Or in other words satisfying the equation above)

For the $\vec{q_{k}}$ we already can simply transform the equations yielded from the matrix multiplication.
The $\vec{r_{k}}$ can be found using the **Gram Schmidt Process**

$$
r_{ij}=\begin{cases}
q_{i}^{*}a_{j} & i<j \\
\lvert \lvert a_{j}-r_{1}q_{1}-\dots-r_{j-1,j}q_{j-1} \rvert  \rvert _{2} & i=j
\end{cases}
$$

# QR Factorization Summary

Given a matrix $A \in \mathbb{R}^{m \times n}$ with linearly independent columns, the **QR factorization** decomposes $A$ into:
$$
A = QR,
$$
where:
- $Q$ is an $m \times n$ matrix with orthonormal columns (obtained via the Gram-Schmidt process).
- $R$ is an $n \times n$ upper triangular matrix with positive diagonal entries.

---

## **Steps to Compute QR Factorization**

### 1. Gram-Schmidt Process:
For each column $\mathbf{a}_i$ of $A$:
1. Compute $\mathbf{v}_i = \mathbf{a}_i - \sum_{j=1}^{i-1} \text{proj}_{\mathbf{u}_j}(\mathbf{a}_i)$.
2. Normalize $\mathbf{v}_i$ to obtain $\mathbf{u}_i = \frac{\mathbf{v}_i}{\|\mathbf{v}_i\|}$.

This gives $Q = [\mathbf{u}_1, \mathbf{u}_2, \dots, \mathbf{u}_n]$.

### 2. Construct $R$:
Two methods are available:
- **Method 1 (Track Operations):** Build $R$ as the sequence of Gram-Schmidt operations "undoes" the transformations.
- **Method 2 (Matrix Multiplication):** Use $R = Q^\top A$.

---

## **Example**

For $A = \begin{bmatrix} -1 & 3 \\ 1 & 5 \end{bmatrix}$:

### Compute $Q$:
1. $\mathbf{v}_1 = \mathbf{a}_1 = \begin{bmatrix} -1 \\ 1 \end{bmatrix}$, $\mathbf{u}_1 = \frac{\mathbf{v}_1}{\|\mathbf{v}_1\|} = \frac{1}{\sqrt{2}} \begin{bmatrix} -1 \\ 1 \end{bmatrix}$.
2. $\mathbf{v}_2 = \mathbf{a}_2 - \text{proj}_{\mathbf{u}_1}(\mathbf{a}_2) = \begin{bmatrix} 3 \\ 5 \end{bmatrix} - \begin{bmatrix} -1 \\ 1 \end{bmatrix} = \begin{bmatrix} 4 \\ 4 \end{bmatrix}$, $\mathbf{u}_2 = \frac{\mathbf{v}_2}{\|\mathbf{v}_2\|} = \frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ 1 \end{bmatrix}$.

Thus:
$$
Q = \frac{1}{\sqrt{2}} \begin{bmatrix} -1 & 1 \\ 1 & 1 \end{bmatrix}.
$$

### Compute $R$:
Using $R = Q^\top A$:
$$
R = \frac{1}{\sqrt{2}} \begin{bmatrix} -1 & 1 \\ 1 & 1 \end{bmatrix} \begin{bmatrix} -1 & 3 \\ 1 & 5 \end{bmatrix} = \sqrt{2} \begin{bmatrix} 1 & 1 \\ 0 & 4 \end{bmatrix}.
$$

---

## **General Formula for $R$**
For $n=3$, $R$ is:
$$
R = \begin{bmatrix}
\|\mathbf{v}_1\| & \langle \mathbf{u}_1, \mathbf{a}_2 \rangle & \langle \mathbf{u}_1, \mathbf{a}_3 \rangle \\
0 & \|\mathbf{v}_2\| & \langle \mathbf{u}_2, \mathbf{a}_3 \rangle \\
0 & 0 & \|\mathbf{v}_3\|
\end{bmatrix}.
$$

Alternatively:
$$
R = Q^\top A.
$$

