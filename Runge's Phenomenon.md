# Runge's Phenomenon

Runge's phenomenon is a problem that occurs when using high-degree polynomials for interpolation over equally spaced points.

## **The Problem**

For a function $f(x)$ defined on an interval $[a, b]$, if we choose $n+1$ equally spaced interpolation points $x_0, x_1, \dots, x_n$, the resulting interpolating polynomial $p_n(x)$ may oscillate significantly, particularly near the endpoints of the interval.

This happens even if $f(x)$ is smooth, such as $f(x) = \frac{1}{1 + x^2}$, the classic Runge example.

## **The Explanation**

1. The oscillations arise because equally spaced points do not minimize the interpolation error.
2. The term $\prod_{k=0}^n (x - x_k)$ in the error bound becomes large near the edges of $[a, b]$.

## **Error Behavior**

If $f(x)$ is smooth and $n \to \infty$, the interpolation error:

$$
\lvert p_n(x) - f(x) \rvert
$$

may grow unbounded, especially near the endpoints.

## **Mitigation**

Runge's phenomenon can be reduced by:

- Using **Chebyshev nodes**, which minimize the maximum value of $\prod_{k=0}^n \lvert x - x_k \rvert$.
- Opting for piecewise polynomial interpolation (splines) instead of a single high-degree polynomial. [[Spline Interpolation]]

## **Key Takeaway**

Runge's phenomenon highlights the limitation of polynomial interpolation with equally spaced points and emphasizes the importance of node distribution.
