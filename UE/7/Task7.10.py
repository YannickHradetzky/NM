import numpy as np

# Define matrix A
A = np.array([
    [1, 2, 3],
    [1, 3, 4], 
    [1, 4, 5],
    [1, 5, 6]
])

# Define vector b
b = np.array([1, 1, 1, 1])

# Compute SVD
U, s, Vh = np.linalg.svd(A, full_matrices=False)

# Compute solution x using SVD
# x = V * Σ⁺ * U^T * b where Σ⁺ is the pseudoinverse of Σ
s_inv = 1/s
x = Vh.T @ (s_inv[:, np.newaxis] * (U.T @ b))

print("Solution x:")
print(x)

# Compute residual norm ||Ax - b||₂
residual = A @ x - b
residual_norm = np.linalg.norm(residual)

print("\nResidual norm ||Ax - b||₂:")
print(residual_norm)

# Add condition number analysis
cond = np.max(s) / np.min(s)
print("\nCondition number:")
print(cond)

# Compute relative error
x_norm = np.linalg.norm(x)
relative_error = residual_norm / np.linalg.norm(b)
print("\nRelative error ||Ax - b||₂ / ||b||₂:")
print(f"{relative_error:.2%}")

# Verify solution quality
print("\nVerification:")
print("||A @ x - b||₂:", residual_norm)
print("||x||₂:", x_norm)
