#NM 
[[Direkt Methods for Systems of Linear equations]]
# LU Factorization
Given a regular matrix A the aim of the LU factorization is to find a lower- and upper triangle matrix $L$ and $U$. Then $Ax=b$ can be solved via first a forward substituition and then a  backward substitution. 
**Gaussian Elimination** 
An algorithm which transforms a given matrix $A$ into an upper-triangular one $U$ by using row operations. For an $n\times n$ matrix it consists of $n-1$ steps. In the $k$ step ($k=1,\dots,n-1$) one eliminates all nonzero entries in the $k$ column below the main diagonal. The connection to the $LU$ factorization becomes apparent when interpreting the row operations as matrix multiplications: Each of the n-1 steps can be realized by multiplication form the left with a lower triangular matrix $L_{k}$.

- $L_k^{-1}$ differs from $L_{k}$ only in the signs of their subdiagonal entries. 
- The product over $k$ simply collects all subdiagonal elements into an identity matrix

$$
\begin{align}
A^{(1)} & =A \\
A^{(k)} & =L_{k-1}A^{(k-1)}=L_{k-1}\dots L_{1}A, \text{ for k=2,...n-1}
\end{align}
$$

So for every $k$ define the numbers (entries of the matrix $L_{k}$)
$$
l_{jk}=\frac{a_{jk}^{(k)}}{a_{kk}^{(k)}}, j=k+1, \dots, n
$$
Then the row operation matrix $L_{k}$ is given by: 
$$
L_k = \begin{bmatrix}
1 &  &  &  &  \\
  & \ddots &  &  &  \\
  &  & 1 &  &  \\
  &  & -\ell_{k+1,k} & 1 &  \\
  &  & \vdots &  & \ddots \\
  &  & -\ell_{n,k} &  &  & 1
\end{bmatrix}
$$
One can define the vector 
$$
l_{k}=(0,\dots,0,l_{k+1,k},\dots,l_{k,n})^{T} \in K^{n}
$$
Then we can use these vectors to rewrite the matrix $L_{k}$:
$$
L_{k}=I-l_{k}e_{k}^{*}
$$
Then the inverse is simply given by reversing the sign of the subdiagonal entries: (Lemma 2.1)
$$
L_{k}^{-1}=I+l_{k}e_{k}^{*}
$$
And the product can be expressed as:  (Lemma 2.2)
$$
L_1^{-1} \cdots L_{n-1}^{-1} = I + \sum_{k=1}^{n-1} (\ell_k e_k^*) =
\begin{bmatrix}
1 &  &  &  \\
\ell_{21} & 1 &  &  \\
\vdots & \ddots & \ddots &  \\
\ell_{n1} & \cdots & \ell_{n,n-1} & 1
\end{bmatrix} = L
$$


## Algorithm
![[Pasted image 20241014111344.png]]


## Theorems

### 2.4)  
[[Skript Numerics 1.pdf]]
To leading order Gaussian elimination requires $\frac{2}{3}n^{3}$ operations.
Solving via [[Forward&Backward Substitution]] requires $\approx2n^{2}$ . 


## Stability
So it is reasonably to us the LU-Factorization as long as for the matrix $A=LU$ if $\kappa(A)\approx \kappa(L)\kappa(U)$


# Pivoting

The stability problem stems from possible divisions with very small numbers in the LU-factorization algorithm. This can be partially circumvented by swapping rows of the matrix $A$. This is fine as long as the same entries are swapped in vector $\vec{b}$. 

in step $k$ in Gaussian elimination we use the matrix element $a_{kk}^{(k)}$ to produce zeros un the $k$-th column.  But we could also use any element $a_{ik}^{(k)}$ , $i>k$, for that. Ideally we want to choose the that element to be the largest in absolut value.  Additionally one can keep the pivot on the main diagonal. This strategy is then called **partial pivoting**. 

So we can write the result of $n-1$ steps of gaussian eliminations as: 
$$
L_{n-1}P_{n-1}\dots L_{1}P_{1}A = U
$$
Where $P_{j}$ are permutation matrices.

**complete pivoting** would be looking for the largest element in all entries $\{ a_{ij}^{k} \}$. This is more expensive and partial pivoting is often sufficient. 

## General Formulas

$A$ is an $n\times n$ matrix. After $n-1$ steps of Gaussian elimination with partial pivoting we can write the resulting upper triangular matrix $U$ in the following way: 

$$
U = L_{n-1}P_{n-1}\dots L_{1}P_{1}A
$$

We can group the matrices in the following way: 
$$
U = L_{-1}PA
$$

**To explain what happened** and understand that we did not just swap matrices in a multiplication:
Lets say : 
$$
\begin{align}
U=L_{3}P_{3}L_{2}P_{2}L_{1}P_{1}A = L_{3}'L_{2}'L_{1}'P_{3}P_{2}P_{1}A
\end{align}
$$
So where $L_k'$ is cleverly chosen to be: 
$$
\begin{align}
L_{3}' & =L_{3} \\
L_{2}' & =P_{3}L_{2}P_{3}^{-1} \\
L_{1}' & =P_{3}P_{2}L_{1}P_{2}^{-1}P_{3}^{-1}
\end{align}
$$

Still we can compute the $L$ matrix of step $k$ in the following way: 

$$
\begin{align}
L_{k}' & =P_{n-1}\dots P_{k+1}\underbrace{ (I-l_{k}e_{k}^{*}) }_{ L_{k}^{-1} }P_{k+1}^{-1}\dots P_{n-1}^{-1} \\
 & =I-\underbrace{ P_{n-1}\dots P_{k+1}l_{k} }_{ l_{k'} }\underbrace{ e_{k}U^{*}P_{k+1}^{-1}\dots P_{n-1}^{-1} }_{ (P_{n-1}\dots P_{k+1}e_{k})^{*} } \\
 & =I-l_{k}'e_{k}^{*}
\end{align}
$$
## Algorithm

![[Pasted image 20241021112537.png]]

## Stability : 

Complicated Topic...see if it is mentioned at all in the lecture 

