import numpy as np
R = np.triu(np.random.rand(50,50)) 
Q = np.linalg.qr(np.random.rand(50,50))[0] 
A = Q@R

Q2, R2 = np.linalg.qr(A)

# Compute forward and backward errors
forward_error_Q = np.linalg.norm(Q - Q2) / np.linalg.norm(Q)
forward_error_R = np.linalg.norm(R - R2) / np.linalg.norm(R)
backward_error = np.linalg.norm(A - Q2@R2) / np.linalg.norm(A)

print(f"Forward error in Q: {forward_error_Q:.2e}")
print(f"Forward error in R: {forward_error_R:.2e}")
print(f"Backward error: {backward_error:.2e}")


