import numpy as np

def back_substitution_upper_triangular(U, b, tol=1e-10, verbose=True):
    """
    Solve an upper triangular system Ux = b using back substitution.
    
    Parameters:
    -----------
    U : numpy.ndarray
        An n x n upper triangular matrix
    b : numpy.ndarray
        A vector of length n containing the right-hand side
    tol : float, optional
        Tolerance for checking diagonal elements (default: 1e-10)
    verbose : bool, optional
        If True, prints verification information (default: True)
        
    Returns:
    --------
    numpy.ndarray
        Solution vector x such that Ux = b
        
    Raises:
    -------
    ValueError
        If U is not square, dimensions don't match, or system is singular
    """
    # Input validation
    if U.shape[0] != U.shape[1]:
        raise ValueError("Matrix U must be square")
    if U.shape[0] != b.shape[0]:
        raise ValueError("Dimensions of U and b must match")
    
    n = U.shape[0]
    x = np.zeros(n)
    
    # Check if the system is singular
    if np.any(np.abs(np.diag(U)) < tol):
        raise ValueError("Matrix U is singular or nearly singular")
    
    # Back substitution
    for i in range(n-1, -1, -1):
        # Calculate sum of known terms
        known_sum = np.dot(U[i, i+1:], x[i+1:])
        # Solve for x[i]
        x[i] = (b[i] - known_sum) / U[i, i]
    
    if verbose:
        _verify_solution(U, x, b)
    
    return x

def _verify_solution(U, x, b, precision=4):
    """
    Helper function to verify and display the solution accuracy.
    
    Parameters:
    -----------
    U : numpy.ndarray
        The upper triangular matrix
    x : numpy.ndarray
        The computed solution
    b : numpy.ndarray
        The original right-hand side
    precision : int, optional
        Number of decimal places to display (default: 4)
    """
    Ux = np.dot(U, x)
    diff = Ux - b
    
    print("\nVerification:")
    print(f"U * x =\n{np.array2string(Ux, precision=precision, separator=', ')}")
    print(f"\nOriginal b =\n{np.array2string(b, precision=precision, separator=', ')}")
    print(f"\nDifference (U*x - b):\n{np.array2string(diff, precision=precision, separator=', ')}")
    print(f"Maximum absolute error: {np.max(np.abs(diff)):.2e}")

# Example usage
if __name__ == "__main__":
    # Example system
    U = np.array([[3, 2, 1],
                  [0, 4, 2],
                  [0, 0, 5]])
    b = np.array([1, 2, 3])
    
    print("Solving the system Ux = b where:")
    print(f"U =\n{U}")
    print(f"b = {b}")
    
    x = back_substitution_upper_triangular(U, b)
    print(f"\nSolution x = {np.array2string(x, precision=4, separator=', ')}")