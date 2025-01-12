import numpy as np

# Generate random 15x30 matrix
np.random.seed(42)  # For reproducibility
A = np.random.randn(15, 30)

# Compute SVD
U, sigma, Vt = np.linalg.svd(A, full_matrices=False)

# Create rank-10 approximation by using only first 10 singular values
sigma_10 = np.copy(sigma)
sigma_10[10:] = 0  # Set smallest 5 singular values to 0

# Construct A10 using truncated singular values
A10 = U @ np.diag(sigma_10) @ Vt

# Calculate Frobenius norm of difference
frob_norm = np.linalg.norm(A - A10, 'fro')

# Calculate sqrt of sum of squares of discarded singular values
sum_squares = np.sqrt(np.sum(sigma[10:]**2))

print(f"Original matrix shape: {A.shape}")
print(f"Singular values: {sigma}")
print(f"\nFrobenius norm ||A - A10||_F: {frob_norm:.6f}")
print(f"sqrt(sum(σ_j^2)) for j=11,...,15: {sum_squares:.6f}")
print(f"\nDifference between measures: {abs(frob_norm - sum_squares):.2e}")
