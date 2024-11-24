import numpy as np
import matplotlib.pyplot as plt

def check_cos_stability(x_range=(-1*np.pi, 1*np.pi), num_points=100000):
    """
    Check the backward stability of numpy's cos implementation by comparing
    the computed result with a perturbed input.
    
    Parameters:
    -----------
    x_range : tuple
        Range of x values to test (default: (-10π, 10π))
    num_points : int
        Number of test points (default: 1000)
        
    Returns:
    --------
    None (displays plot and prints statistics)
    """
    # Generate test points
    x = np.linspace(x_range[0], x_range[1], num_points)
    
    # Compute cos(x) using numpy
    y = np.cos(x)
    
    # Machine epsilon
    eps = np.finfo(float).eps
    
    # Introduce small perturbation to input
    delta_x = eps * x
    x_perturbed = x + delta_x
    y_perturbed = np.cos(x_perturbed)
    
    # Compute relative errors
    rel_input_error = np.abs(delta_x / x)
    rel_output_error = np.abs((y_perturbed - y) / y)
    
    # Compute condition number |x*sin(x)/cos(x)|
    condition_number = np.abs(x * np.sin(x) / np.cos(x))
    
    # Plot results
    plt.figure(figsize=(12, 8))
    plt.semilogy(x, rel_input_error, label='Relative input perturbation')
    plt.semilogy(x, rel_output_error, label='Relative output error')
    plt.semilogy(x, condition_number * eps, label='Expected error bound')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('Relative error')
    plt.title('Backward stability analysis of np.cos')
    plt.legend()
    
    # Print statistics
    print(f"Maximum relative input perturbation: {np.max(rel_input_error):.2e}")
    print(f"Maximum relative output error: {np.max(rel_output_error):.2e}")
    print(f"Maximum condition number: {np.max(condition_number):.2e}")
    
    plt.show()

# Run the stability check
if __name__ == "__main__":
    check_cos_stability()
