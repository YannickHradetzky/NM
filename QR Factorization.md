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

