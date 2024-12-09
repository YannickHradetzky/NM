# Task 1
```python 
import numpy as np

def fourier_matrix(n):
    """
    Generate the Fourier matrix F_n of size n x n.
    
    Args:
    - n: Size of the Fourier matrix
    
    Returns:
    - F_n: Fourier matrix
    """
    omega = np.exp(-2j * np.pi / n)  # nth root of unity
    indices = np.arange(n)
    F_n = np.power(omega, np.outer(indices, indices)) / np.sqrt(n)  # Normalize
    return F_n

def compute_condition_numbers(F_n):
    """
    Compute the condition numbers of the Fourier matrix F_n.
    
    Args:
    - F_n: Fourier matrix
    
    Returns:
    - kappa_1: Condition number in L1 norm
    - kappa_2: Condition number in L2 norm
    - kappa_inf: Condition number in Linf norm
    """
    F_n_inv = np.linalg.inv(F_n)
    norm_1 = np.linalg.norm(F_n, 1) * np.linalg.norm(F_n_inv, 1)
    norm_2 = np.linalg.norm(F_n, 2) * np.linalg.norm(F_n_inv, 2)
    norm_inf = np.linalg.norm(F_n, np.inf) * np.linalg.norm(F_n_inv, np.inf)
    return norm_1, norm_2, norm_inf

def recursive_factorization(n):
    """
    Recursive factorization of Fourier matrix F_n for n even.
    
    Args:
    - n: Size of the Fourier matrix (must be even)
    
    Returns:
    - F_n_recursive: Recursively factorized Fourier matrix
    """
    if n % 2 != 0:
        raise ValueError("n must be even for recursive factorization")
    
    # Base case: 2x2 Fourier matrix
    if n == 2:
        return np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    
    # Recursive step
    m = n // 2
    F_m = recursive_factorization(m)
    omega = np.exp(-2j * np.pi / n)
    Omega = np.diag([omega**k for k in range(m)])
    
    # Block structure
    top = np.hstack((F_m, Omega @ F_m))
    bottom = np.hstack((F_m, -Omega @ F_m))
    return np.vstack((top, bottom))

def unitary_factorization(F_n):
    """
    Compute the unitary factorization of the Fourier matrix F_n.
    
    Args:
    - F_n: Fourier matrix
    
    Returns:
    - Q: Unitary matrix
    - D: Diagonal matrix
    """
    eigenvalues, Q = np.linalg.eig(F_n)
    D = np.diag(eigenvalues)
    return Q, D

# Solve for an example size n
n = 4
F_n = fourier_matrix(n)

# Compute condition numbers
kappa_1, kappa_2, kappa_inf = compute_condition_numbers(F_n)

# Recursive factorization
F_n_recursive = recursive_factorization(n)

# Unitary factorization
Q, D = unitary_factorization(F_n)

# Display results
print(f"Fourier matrix F_{n}:\n{F_n}\n")
print(f"Condition Numbers:\nL1 norm: {kappa_1}\nL2 norm: {kappa_2}\nLinf norm: {kappa_inf}\n")
print(f"Recursive Factorization of F_{n}:\n{F_n_recursive}\n")
print(f"Unitary Factorization of F_{n}:\nQ (Unitary Matrix):\n{Q}\nD (Diagonal Matrix):\n{D}\n")
```


# Task 2
```python 
import numpy as np
import matplotlib.pyplot as plt

def compute_dft_coefficients(n, y):
    # Compute x_j
    x_j = 2 * np.pi * np.arange(n) / n
    
    # Perform normalized FFT
    y_fft = np.fft.fft(y) / n
    a_0 = y_fft[0].real  # First coefficient
    
    # Extract coefficients
    a_k = 2 * y_fft[1:(n//2)].real  # Cosine coefficients
    b_k = -2 * y_fft[1:(n//2)].imag  # Sine coefficients
    
    if n % 2 == 0:
        a_k = np.append(a_k, y_fft[n//2].real)  # Append a_n/2 for even n
    
    # Define interpolant function
    def interpolant(x):
        result = a_0
        for k in range(1, len(a_k) + 1):
            result += a_k[k-1] * np.cos(k * x)
            if k <= len(b_k):
                result += b_k[k-1] * np.sin(k * x)
        return result
    
    return x_j, a_0, a_k, b_k, interpolant

def plot_interpolant(x_j, y, interpolant, title):
    x_plot = np.linspace(0, 2 * np.pi, 500)
    y_plot = interpolant(x_plot)
    
    plt.figure(figsize=(8, 6))
    plt.plot(x_plot, y_plot, label='Interpolant', color='blue')
    plt.scatter(x_j, y, color='red', label='Data Points', zorder=5)
    plt.title(title)
    plt.xlabel('$x$')
    plt.ylabel('$P(x)$')
    plt.legend()
    plt.grid()
    plt.show()

# Data for the given cases
cases = {
    "n=4, y=(1, 0, 1, 0)": (4, np.array([1, 0, 1, 0])),
    "n=4, y=(0, 1, 0, 1)": (4, np.array([0, 1, 0, 1])),
    "n=8, y=(1, 0, 1, 0, 1, 0, 1, 0)": (8, np.array([1, 0, 1, 0, 1, 0, 1, 0])),
    "n=8, y=(0, 1, 0, 1, 0, 1, 0, 1)": (8, np.array([0, 1, 0, 1, 0, 1, 0, 1]))
}

# Solve and display results for each case
for case_name, (n, y) in cases.items():
    print(f"Case: {case_name}")
    
    # Compute coefficients and interpolant
    x_j, a_0, a_k, b_k, interpolant = compute_dft_coefficients(n, y)
    
    # Display coefficients
    print(f"x_j: {x_j}")
    print(f"a_0: {a_0}")
    print(f"a_k: {a_k}")
    print(f"b_k: {b_k}")
    print()
    
    # Plot interpolant
    plot_interpolant(x_j, y, interpolant, title=f"Trigonometric Interpolant for {case_name}")
```