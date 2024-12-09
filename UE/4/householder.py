import numpy as np

def solve_qr_householder(A, b):
    Q,R = np.linalg.qr(A)
    # if m > n, compute least squares solution
    print(f"Condition number of A: {np.linalg.cond(A)}")
    if A.shape[0] == A.shape[1]:
        print("System is square")
        # solve Ax = b using QR factorization
        x_with_A = np.linalg.solve(A, b)
        # instead of Ax = b we use QR factorization so we solve Rx = Q^T @ b
        x_with_Q = np.linalg.solve(R, Q.T @ b)
    else:
        print("System is not square")
        print("Using least squares solution")
        # solve least squares problem
        x_with_A = np.linalg.lstsq(A, b, rcond=None)[0]
        # instead of Ax = b we use QR factorization so we solve Rx = Q^T @ b
        x_with_Q = np.linalg.lstsq(R, Q.T @ b, rcond=None)[0]

    # print both solutions and the relative and absolute errors
    print(f"x_with_A: {x_with_A}")
    print(f"x_with_Q: {x_with_Q}")
    print(f"Relative error: {np.linalg.norm(x_with_A - x_with_Q) / np.linalg.norm(x_with_A)}")
    print(f"Absolute error: {np.linalg.norm(x_with_A - x_with_Q)}")
    return x_with_A, x_with_Q

A = np.random.rand(4,3)
b = np.random.rand(4)
x_with_A, x_with_Q = solve_qr_householder(A, b)