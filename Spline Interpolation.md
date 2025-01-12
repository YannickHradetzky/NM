[[Skript Numerics 1.pdf]]
# Spline Interpolation

Spline interpolation is a method of fitting a piecewise polynomial function to a set of data points. Unlike high-degree polynomial interpolation, splines provide a smooth approximation without the oscillatory behavior seen in [[Runge's Phenomenon]].

## **The Idea**

Given $n+1$ data points $(x_0, y_0), \dots, (x_n, y_n)$, spline interpolation constructs a piecewise polynomial $s(x)$ such that:

1. $s(x_i) = y_i$ for all $i \in \{0, \dots, n\}$.
2. $s(x)$ is smooth, typically ensuring continuity of the first and second derivatives at the boundaries of the pieces.

## **Types of Splines**

- [[Linear Splines]]: Piecewise linear functions with continuity at the data points.
- **Cubic Splines**: Piecewise cubic polynomials with continuity of the first and second derivatives.
- **Higher-Order Splines**: Less commonly used due to increased complexity.

## **Cubic Splines**

Cubic splines are the most commonly used. For each interval $[x_i, x_{i+1}]$, the spline is a cubic polynomial:

$$
s(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
$$

The coefficients $a_i, b_i, c_i, d_i$ are determined by:

1. The interpolation condition: $s(x_i) = y_i$ and $s(x_{i+1}) = y_{i+1}$.
2. Smoothness conditions: $s'(x)$ and $s''(x)$ are continuous.
3. Boundary conditions, such as natural spline ($s''(x_0) = s''(x_n) = 0$) or clamped spline (specified derivatives at endpoints).

## **Advantages**

- Avoids the oscillations of high-degree polynomial interpolation.
- Provides a smooth approximation.
- Local changes in the data affect only nearby segments of the spline.

## **Applications**

Spline interpolation is widely used in:

- Computer graphics for curve modeling.
- Data fitting in engineering and scientific computations.
- Motion planning in robotics.

## **Key Takeaway**

Spline interpolation divides the interval into smaller subintervals and fits low-degree polynomials to each, ensuring smoothness and stability across the entire range.

# General Stuff from Script

## Lemma 3.2 

Let $n,p \in N$ and $\Delta$ a partition of $[a,b]$.

1) The Set of Splines of order $p$ $S^{P}(\Delta)$ is a vector space of dimension $n+p$
2) Let $0\leq k\leq p-1$. If $s \in S^{p}(\Delta)$, then $s ^{(k)}\in S^{p-k}$. So the $k$-th derivative of a $p-th$ order spline is a spline of order $p-k$.

