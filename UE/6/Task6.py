import numpy as np

def composite_trapezoidal(x, y):
    """
    Implements the composite trapezoidal rule for numerical integration.
    
    Parameters:
    x (array): Vector of nodes in ascending order (x[i] < x[i+1])
    y (array): Vector of function values corresponding to x
    
    Returns:
    float: Estimate of the definite integral
    """
    # Input validation
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x) < 2:
        raise ValueError("At least 2 points are required")
    if not np.all(np.diff(x) > 0):
        raise ValueError("x must be strictly increasing")
    
    # Calculate step sizes between consecutive x values
    dx = np.diff(x)
    
    # Apply trapezoidal rule to each subinterval and sum
    # For each subinterval [x[i], x[i+1]]:
    # Area ≈ (y[i] + y[i+1])/2 * (x[i+1] - x[i])
    integral = np.sum(0.5 * (y[:-1] + y[1:]) * dx)
    
    return integral

# Example usage
if __name__ == "__main__":
    # Test with a simple function, e.g., f(x) = x^2 from 0 to 1
    x = np.linspace(0, 1, 11)  # 11 points including endpoints
    y = x**2  # function values
    
    result = composite_trapezoidal(x, y)
    exact = 1/3  # exact value of integral of x^2 from 0 to 1
    
    print(f"Numerical result: {result:.6f}")
    print(f"Exact result: {exact:.6f}")
    print(f"Absolute error: {abs(result - exact):.6e}")
