import numpy as np
import matplotlib.pyplot as plt
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
    integral = 0
    for i in range(len(x) - 1):
        integral += 0.5 * (y[i] + y[i+1]) * dx[i]
    return integral



# Example usage
if __name__ == "__main__":
    # Test with a simple function, e.g., f(x) = x^2 from 0 to 1
    x = np.linspace(0, 1, 11)  # 11 points including endpoints
    y = x**2  # function values

    # plot the function and the trapezoidal rule
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='f(x) = x^2')
    plt.plot(x, y, 'ro', label='Data points')
    
    # Draw the trapezoidal blocks
    for i in range(len(x)-1):
        # Create vertices of the trapezoid
        trap_x = [x[i], x[i], x[i+1], x[i+1]]
        trap_y = [0, y[i], y[i+1], 0]
        plt.fill(trap_x, trap_y, alpha=0.2)  # alpha controls transparency
        
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Trapezoidal Rule')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    result = composite_trapezoidal(x, y)
    print(result)
    
    print(f"Numerical result: {result:.6f}")


