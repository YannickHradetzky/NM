#NM 

# 1) Linear Independency and Gram-Schmidt process

## Linearly independent?
To check if the columns of the matrix are linear independent we can compute the determinant of the matrix is non zero:

$$
\begin{align}
\det[\begin{bmatrix}
1 & 2 & 3 \\
0 & 2 & 3 \\
0 & 0 & 3
\end{bmatrix}]  & = [(1*2*3) + (2*3*0)+(0*0*3)]-[(0*2*3)+(0*2*3)+(0*3*1)] \\
 & =[6+0+0]-[0] = 6 \neq 0
\end{align}
$$
So we have verified that the columns of the matrix are indeed linearly independent.

## Gram-Schmidt process: 

Used to obtain orthonormal vectors that are easier to work with but span the same space as the original column vectors

Take the first vector $\vec{v_{1}}=(1,0,0)^{T}$ and normalize it. But it is already normalized so we gonna refer to it as $\vec{u_{1}}$

Now we take the second vector $\vec{v_{2}} = (2,2,0)^{T}$ and compute its projection onto $\vec{u_{1}}$.

$$
\begin{align} 
\text{proj}_{\vec{u_{1}}}(\vec{v_{2}})  & = \frac{\vec{v_{2}}\vec{\cdot u_{1}}}{||\vec{u_{1}}||_{2}}\vec{\cdot u_{1}} \\
	 & = \frac{1*2+2*0+0*0}{1}\cdot \begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix} = \begin{bmatrix}
2 \\
0 \\
0
\end{bmatrix}

\end{align}
$$

Next we take this projection and subtract it from $\vec{v_{2}}$

$$
\vec{v_{2}}' = \vec{v_{2}}- \text{proj}_{\vec{u_{1}}}(\vec{v_{2}}) =  \begin{bmatrix}
2 \\
2 \\
0 
\end{bmatrix} - 
\begin{bmatrix}
2 \\
0 \\
0
\end{bmatrix} = 
\begin{bmatrix}
0 \\
2 \\
0
\end{bmatrix}
$$
Then we normalize it by dividing by its $||\dots||_{2}$ Norm which is clearly $2$ and we yield $\vec{u_{2}}=(0,1,0)^{T}$

For the third vector $\vec{v_{3}}=(3,3,3)^{T}$ we calculate its projection onto $\vec{u_{1}}$ and $\vec{u_{2}}$ and subtract both from $\vec{v_{3}}$:
$$
\begin{align}
\text{proj}_{\vec{u_{1}}}(\vec{v_{3}})  & = \frac{\vec{v_{3}}\vec{\cdot u_{1}}}{||\vec{u_{1}}||_{2}}\vec{\cdot u_{1}} \\
 & = 3 \cdot \begin{bmatrix}
1 \\
0 \\
0 \\
\end{bmatrix} = \begin{bmatrix}
3 \\
0 \\
0
\end{bmatrix}
\end{align}
$$

and 
$$
\begin{align}
\text{proj}_{\vec{u_{2}}}(\vec{v_{3}})  & = \frac{\vec{v_{3}}\vec{\cdot u_{2}}}{||\vec{u_{2}}||_{2}}\vec{\cdot u_{2}} \\
 & = 3 \cdot \begin{bmatrix}
0 \\
1 \\
0 \\
\end{bmatrix} = \begin{bmatrix}
0 \\
3 \\
0
\end{bmatrix}
\end{align}
$$

We subtract it from $\vec{v_{3}}$:
$$
\vec{v_{3}}' = \vec{v_{3}} - \text{proj}_{\vec{u_{2}}}(\vec{v_{3}}) - \text{proj}_{\vec{u_{1}}}(\vec{v_{3}}) = \begin{bmatrix}
0 \\
0 \\
3
\end{bmatrix}
$$

We normalize it by dividing by its norm ($3$) and yield the last orthonormal vector $\vec{u_{3}}=(0,0,1)^{T}$

So in total we get the set of vectors: 

$$
\{ \vec{u_{1}}, \vec{u_{2}}, \vec{u_{3}} \} = \{ \begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}, \begin{bmatrix}
0 \\
1 \\
0
\end{bmatrix}, \begin{bmatrix}
0 \\
0 \\
1
\end{bmatrix} \}
$$
Which are clearly orthonormal and process the same Span as the original vectors but are easier to work with

# 2) reflection matrix 

Definition of the reflection matrix : 
$$
H = I - \frac{2}{w^*w}ww^{*}
$$
## a) H is self-adjoint (hermitian) and unitary/orthogonal and involutory

**Hermitian**
For that we compute $H^{*}$

$$
\begin{align}
H^{*}  & = \left( I - \frac{2}{\underbrace{ w^*w }_{ \in R }}\underbrace{ ww^{*} }_{ \in K^{n\times n} } \right)^* \\
 & = I - \frac{2}{w^{*}w}(ww^{*})^{*} =  I - \frac{2}{w^{*}w}(w^{*})^{*}w^{*} \\
 & = I - \frac{2}{w^{*}w}ww^{*} = H
\end{align}
$$
**unitary/orthogonal** and **involutory**
For that we have to show $H^{*}H = I$. Since we already know that $H$ is hermitian this simplifies to $H^{2}$.

$$
\begin{align}
H^{2} & =\left( I - \frac{2}{w^*w}ww^{*} \right)\cdot\left( I - \frac{2}{w^*w}ww^{*} \right) \\
 & = I^{2} - 2\cdot I\cdot\frac{2}{w^*w}ww^{*} + \left( \frac{2}{w^*w}ww^{*} \right)^{2} \\
 & =I^{2}-4\frac{ww^{*}}{w^{*}w}+ 4 \underbrace{ \left( \frac{ww^{*}}{w^{*}w}  \right)^{2} }_{ \frac{w(w^{*}w)w^{*}}{(w^{*}w)^{2}} = \frac{ww^{*}}{w^{*}w} } \\
 & =I^{2}=I
\end{align}
$$

## b) Determine eigenvalues, determinant and singular values

Luckily there is no need for computing the eigenvalues. Since the reflection matrix reflects across the subspace orthogonal to the vector $w$. So any vector in the direction of $w$ will be "flipped" with a eigenvalue of $-1$. Any vector in a orthogonal subspace to $w$ will remain unchanged with an eigenvalue $+1$. 

$$
\lambda_1 = -1 \quad \text{(with multiplicity 1)}
$$
$$

\lambda_2 = 1 \quad \text{(with multiplicity  n-1 )},

$$

The determinant can be computed by calculating the product of the eigenvalues. 
$$
\det(H) = \prod_{i}^{n}\lambda_{i} = -1 \cdot 1^{(n-1)} = -1
$$

The singular values of a matrix  $H$  are the square roots of the eigenvalues of the matrix  $H^* H$ . Since we already showed that  $H$  is **unitary** (or **orthogonal**), we have  $H^* H = I$ . This means all eigenvalues of  $H^* H$  are  $1$ , and thus the singular values of  $H$ are all  $1$ .

Therefore, the **singular values of  H**  are:
$$
\sigma_1 = \sigma_2 = \cdots = \sigma_n = 1
$$

# 3)  Norm (in)equalities
##  **$||x||_{\infty}\leq||x||_{2}$

The infinity norm picks the largest element of $x \in K^{n}$ while the Euclidean norm  $|x\|2 = \left(\sum{i=1}^{n} |x_i|^2\right)^{1/2}$
Therefore the $||x||_{2}$ is always as least as large as the $||x||_{\infty}$ norm. 
Equality holds if all the elements have the same value e.g: $(1,0,0,\dots,0)^{T}$

## $||x||_{2}\leq \sqrt{ n }||x||_{\infty}$

This holds because the $||x||_{\infty}$ norm gives us an upper bound on the elements that get squared and summed. So for sure $n||x||_{\infty}\geq \sum_{i=1}^{n}x_{i}$. If we then look at the sum of squares: 

$$
\begin{align}
n||x||_{\infty}  & \geq \sum_{i=1}^{n}x_{i}^{2} \\
\sqrt{ n }||x||_{\infty}  & \geq \sqrt{ \sum_{i=1}^{n}x_{i}^{2} }
\end{align}
$$

Equality will hold if all elements are the same e.g: $(1,1,\dots,1)^{T}$

# 4) Norm limits

The Definition of the $p$-Norm is: 

$$
||x||_{p}= \left( \sum_{i}^{n}|x_{i}|^{p} \right)^{1/p}
$$
## Intuition 
The larger $p$ gets the more the sum will get dominated by the $p$-th power of the largest term which is by definition the $||x||_{\infty}$. Lets assume $||x||_{\infty}=x_{k}$

So for very large $p$:
$$
\begin{align}
||x||_{p}= \left( \sum_{i}^{n}|x_{i}|^{p} \right)^{1/p}\approx x_{k}
\end{align}
$$

And as $p$ gets bigger and bigger the domination of the largest term will become stronger and stronger. Therefore taking the limit does indeed converge to the largest element or in other words the $\infty$-Norm. 

# 5) Frobenius Norm $||A||_{F}$

We have to show that:

$$
||A||_{F} = (\mathrm{Tr}(A^{*}A))^{1/2} = (\mathrm{Tr}(AA^{*}))^{1/2}
$$

With $||A||_{F}$ being defined as follows:
$$
||A||_{F}=\sqrt{ \sum_{i}^{d}\sum_{k}^{d}|A_{ik}|^{2} }
$$

In general the trace has the following property : $\mathrm{Tr}(AB)=\mathrm{Tr}(BA)$ and is defined as the sum over the diagonal entries $\mathrm{Tr}A=\sum_{i}A_{ii}$

Now lets look at $A^{*}A$ more closely: 
The conjugate transpose works as follows: $(A_{ij})^{*}=\bar{A_{ji}}$
$$
(A^{*}A)_{ij} = \sum_{k=1}^{d}(A_{ik})^{*}A_{kj} = \sum_{k=1}^{d}A_{ki}A_{kj}
$$

The diagonal entries are $(A^{*}A)_{ii}=\sum_{k=1}^{d}(A_{ik})^{*}A_{ki}=\sum_{k=1}^{d}|A_{ki}|^{2}$

Putting this all together the trace becomes: 

$$
\begin{align}
\mathrm{Tr}(A^{*}A) & = \sum_{i=1}^{d}(A^{*}A)_{ii} \\
 & = \sum_{i=1}^{d}\sum_{k=1}^{d}|A_{ik}|^{2}
\end{align}
$$
So we have shown that $||A||_{F}^{2}=\sum_{i=1}^{d}\sum_{k=1}^{d}|A_{ik}|^{2}$. A trivial application of the sqare root will yield the desired result. 
In principle the $\mathrm{Tr}$ operation is invariant under cyclic permutations. But the proof for the second equality would follow the same outline as seen above 

# 6) Condition numbers $\kappa_{1},\kappa_{2}$ and $\kappa_{\infty}$

The **condition number** of a matrix measures how sensitive the solution of a linear system is to small changes or errors in the input. A high condition number indicates that even small errors in the input (such as round-off errors) could lead to large errors in the solution. The condition number depends on the matrix norm being used.

The definition of the $p$-kondition number is: 
$$
\kappa_{p}=||A||_{p}||\cdot||A^{-1}||_{p}
$$
So lets first compute the inverse of the matrix. A quick check shows that the determinant of the matrix is not zero.

$$
A = \begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix}
$$
$$
A^{-1} = \frac{1}{\det(A)}\begin{bmatrix}
2 & 1 \\
1 & 2
\end{bmatrix} = \frac{1}{3}\begin{bmatrix}
2 & 1 \\
1 & 2
\end{bmatrix}
$$
## $p=1$

$$
||A||_{1} = \max_{1\leq j\leq n}\sum_{i=1}^{m}|A_{ij}|
$$
which is simply the maximum absolute column sum of the matrix. 

$$
\kappa_{1}=||A||_{1}||\cdot||A^{-1}||_{1} = 3 \cdot 1 = 3
$$

## $p=2$
$$
\lvert \lvert A \rvert  \rvert _{2} = \sqrt{ \lambda_{max}(A^{*}A) }=\sigma_{max}(A)
$$

**The singular values of a matrix  $H$  are the square roots of the eigenvalues of the matrix  $H^* H$ .

So lets calculate $A^{*}$. Which is not difficult since all values are real. So all we have to do is to transpose the matrix.

$$
A^{*}=
\begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix} = A
$$
Which means that the matrix is hermitian. So $A^{*}A=A^{2}$
$$
A^{2}=\begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix}\cdot \begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix}=\begin{bmatrix}
5 & -4 \\
-4 & 5
\end{bmatrix}
$$
Let's calculate the eigenvalues by solving $\det(A^{2}-I\lambda)$

$$
\begin{align}
\det(A^{2}-I\lambda)=\det\begin{bmatrix}
5-\lambda & -4 \\
-4 & 5-\lambda
\end{bmatrix}  & = (5-\lambda)^{2} - 16 =0\\
 & =5-\lambda = \pm 4
\end{align}
$$
So we find two solutions

$$
\lambda=\begin{cases}
\lambda_{1} = 1 \\
\lambda_{2} = 9
\end{cases}
$$

So we find that the $\lvert \lvert A \rvert \rvert_{2}=\sqrt{ 9 }=3$

## $p=\infty$

$$
\lvert \lvert A \rvert  \rvert _{\infty}=\max_{1\leq i\leq m}\sum_{i=1}^{n}\lvert A_{ij} \rvert 
$$
which is simply the maximum absolute row sum of the matrix.

So we get: 
$$
\begin{align}
\kappa_{\infty} & =\lvert \lvert A \rvert  \rvert _{\infty}\cdot \lvert \lvert A^{-1} \rvert  \rvert _{\infty} \\
 & = 3 \cdot 1 = 3
\end{align}
$$

# 7) Conditioning and perturbation

We solve the unperturbed system $Ax=b$

**unperturbed system**
First we bring the matrix in upper triangle form: 

$$
\begin{align}
\begin{bmatrix}
1 & 4 \\
1 & 4.001
\end{bmatrix} \cdot
\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} = \begin{bmatrix}
4 \\
4
\end{bmatrix} \\
\begin{bmatrix}
1 & 4 \\
0 & 0.001
\end{bmatrix}\cdot\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} = \begin{bmatrix}
4 \\
0
\end{bmatrix}
\end{align}
$$
Which yields the following two equations:
$$
\begin{align}
x_{1} + 4x_{2}=4 \\
0.001x_{2} = 0
\end{align}
$$
So we have the solution: 
$$
x = \begin{bmatrix}
4 \\
0
\end{bmatrix}
$$
**perturbed system**

$$
\begin{align}
\begin{bmatrix}
1 & 4 \\
1 & 4.001
\end{bmatrix} \cdot
\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} = \begin{bmatrix}
4 \\
4
\end{bmatrix} + \begin{bmatrix}
0 \\
0.001
\end{bmatrix} \\
\begin{bmatrix}
1 & 4 \\
0 & 0.001
\end{bmatrix}\cdot\begin{bmatrix}
x_{1} \\
x_{2}
\end{bmatrix} = \begin{bmatrix}
4 \\
0.001
\end{bmatrix}
\end{align}
$$
So we have the following system of equations:

$$
\begin{align}
x_{1} + 4x_{2}=4 \\
0.001x_{2} = 0.001
\end{align}
$$
Which yields the solution 
$$
x'=\begin{bmatrix}
0 \\
1
\end{bmatrix}
$$
**relative error**
There is some theorem that all norms are equivalent so we are gonna use the $\lvert \lvert \dots \rvert \rvert_{2}$ norm
$$
\begin{align}
\mathcal{E}_{rel} & = \frac{\lvert \lvert x'-x \rvert  \rvert_{2}}{\lvert \lvert x \rvert  \rvert_{2}} \\
 & =\frac{\sqrt{ (-4)^{2}+1^{2} }}{\sqrt{ 16 }}=\frac{\sqrt{ 17 }}{4} \approx 1.03
\end{align}
$$
This means that the change in the solution is of the same magnitude as the original solution. Therefore the problem/linear-system is ill-conditioned. 



