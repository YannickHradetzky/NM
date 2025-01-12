import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

def generate_symmetric_matrix(n):
    """Generate a random symmetric matrix."""
    A = np.random.randn(n, n)
    return (A + A.T) / 2

def power_iteration(A, max_iter=100, tol=1e-10):
    """
    Implement power iteration method with Rayleigh quotient.
    Returns convergence history of eigenvalue and eigenvector errors.
    """
    n = A.shape[0]
    
    # Get true eigenvalues and eigenvectors for error calculation
    eigenvalues, eigenvectors = eig(A)
    idx = np.argsort(abs(eigenvalues))[::-1]
    true_eigenvalues = eigenvalues[idx]
    true_eigenvector = eigenvectors[:, idx[0]]
    
    # Initialize random vector
    x = np.random.randn(n)
    x = x / np.linalg.norm(x)
    
    # Initialize storage for convergence history
    eigenvalue_errors = []
    eigenvector_errors = []
    theoretical_errors = []
    
    # Calculate convergence rate q = |λ₂/λ₁|
    q = abs(true_eigenvalues[1] / true_eigenvalues[0])
    
    for k in range(max_iter):
        # Power iteration step
        y = A @ x
        
        # Rayleigh quotient
        rayleigh = (x.T @ A @ x) / (x.T @ x)
        
        # Normalize
        x_new = y / np.linalg.norm(y)
        
        # Ensure consistent direction with true eigenvector
        if np.dot(x_new, true_eigenvector) < 0:
            x_new = -x_new
        
        # Calculate errors
        eigenvalue_error = abs(rayleigh - true_eigenvalues[0])
        eigenvector_error = np.linalg.norm(x_new - true_eigenvector)
        theoretical_error = q ** k
        
        # Store errors
        eigenvalue_errors.append(eigenvalue_error)
        eigenvector_errors.append(eigenvector_error)
        theoretical_errors.append(theoretical_error)
        
        # Update vector
        x = x_new
        
        # Check convergence
        if eigenvalue_error < tol:
            break
    
    return eigenvalue_errors, eigenvector_errors, theoretical_errors, q

def plot_convergence(eigenvalue_errors, eigenvector_errors, theoretical_errors, q):
    """Plot convergence behavior in semi-logarithmic scale."""
    iterations = range(len(eigenvalue_errors))
    
    plt.figure(figsize=(10, 6))
    plt.semilogy(iterations, eigenvalue_errors, 'b-', label='Eigenvalue Error')
    plt.semilogy(iterations, eigenvector_errors, 'r-', label='Eigenvector Error')
    plt.semilogy(iterations, theoretical_errors, 'g--', label=f'Theoretical (q^k), q={q:.4f}')
    
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Power Iteration Convergence')
    plt.legend()
    plt.tight_layout()
    
    return plt

def main():
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate symmetric matrix
    n = 5
    A = generate_symmetric_matrix(n)
    
    # Run power iteration
    eigenvalue_errors, eigenvector_errors, theoretical_errors, q = power_iteration(A)
    
    # Plot results
    plt = plot_convergence(eigenvalue_errors, eigenvector_errors, theoretical_errors, q)
    
    # Print final errors
    print(f"Final eigenvalue error: {eigenvalue_errors[-1]:.2e}")
    print(f"Final eigenvector error: {eigenvector_errors[-1]:.2e}")
    print(f"Convergence rate (q = |λ₂/λ₁|): {q:.4f}")
    
    plt.show()

if __name__ == "__main__":
    main()