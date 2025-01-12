import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig, solve

def inverse_iteration(A, mu, max_iter=100, tol=1e-10, use_rayleigh=False):
    """
    Implement inverse iteration method with optional Rayleigh quotient shift.
    Returns convergence history of eigenvalue and eigenvector errors.
    """
    n = A.shape[0]
    
    # Get true eigenvalues and eigenvectors for error calculation
    eigenvalues, eigenvectors = eig(A)
    
    # Initialize random vector
    x = np.random.randn(n)
    x = x / np.linalg.norm(x)
    
    # Initialize storage for convergence history
    eigenvalue_errors = []
    eigenvector_errors = []
    
    # Find closest eigenvalue/vector to shift
    idx = np.argmin(abs(eigenvalues - mu))
    true_eigenvalue = eigenvalues[idx]
    true_eigenvector = eigenvectors[:, idx]
    
    # Ensure consistent direction
    if np.dot(x, true_eigenvector) < 0:
        true_eigenvector = -true_eigenvector
    
    shift = mu
    for k in range(max_iter):
        # Solve linear system
        y = solve(A - shift * np.eye(n), x)
        
        # Normalize
        x_new = y / np.linalg.norm(y)
        
        # Calculate Rayleigh quotient
        rayleigh = (x_new.T @ A @ x_new) / (x_new.T @ x_new)
        
        # Update shift if using Rayleigh quotient
        if use_rayleigh:
            shift = rayleigh
            
        # Ensure consistent direction
        if np.dot(x_new, true_eigenvector) < 0:
            x_new = -x_new
            
        # Calculate errors
        eigenvalue_error = abs(rayleigh - true_eigenvalue)
        eigenvector_error = np.linalg.norm(x_new - true_eigenvector)
        
        # Store errors
        eigenvalue_errors.append(eigenvalue_error)
        eigenvector_errors.append(eigenvector_error)
        
        # Update vector
        x = x_new
        
        # Check convergence
        if eigenvalue_error < tol:
            break
            
    return eigenvalue_errors, eigenvector_errors, rayleigh, x

def plot_convergence(errors, title):
    """Plot convergence behavior in semi-logarithmic scale."""
    iterations = range(len(errors))
    
    plt.semilogy(iterations, errors, 'b-', label='Error')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title(title)
    plt.legend()

def main():
    # Define test matrix
    A = np.array([[3, 2, 1],
                  [2, 3, 1],
                  [1, 2, 3]])
    
    # Test shifts
    shifts = [5, 7]
    
    plt.figure(figsize=(15, 10))
    
    # Test with standard inverse iteration
    for i, mu in enumerate(shifts, 1):
        eigenvalue_errors, eigenvector_errors, final_eigenvalue, final_eigenvector = inverse_iteration(A, mu)
        
        plt.subplot(2, 2, i)
        plot_convergence(eigenvalue_errors, f'Standard Inverse Iteration (μ={mu})')
        
        print(f"\nResults for shift μ={mu}:")
        print(f"Final eigenvalue: {final_eigenvalue:.6f}")
        print(f"Final eigenvector: {final_eigenvector}")
        
    # Test with Rayleigh quotient shift
    for i, mu in enumerate(shifts, 3):
        eigenvalue_errors, eigenvector_errors, final_eigenvalue, final_eigenvector = inverse_iteration(A, mu, use_rayleigh=True)
        
        plt.subplot(2, 2, i)
        plot_convergence(eigenvalue_errors, f'Rayleigh Quotient Shift (initial μ={mu})')
        
        print(f"\nResults for Rayleigh quotient shift (initial μ={mu}):")
        print(f"Final eigenvalue: {final_eigenvalue:.6f}")
        print(f"Final eigenvector: {final_eigenvector}")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
