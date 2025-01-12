import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Read and convert image to grayscale numpy array
img = Image.open('/Users/yhra/Documents/Master/Semester_3/NM/UE/7/IMG_6049.jpeg').convert('L')
A = np.array(img)

# Perform SVD
U, sigma, Vt = np.linalg.svd(A, full_matrices=False)

# Calculate number of singular values for each percentage
n_singular_values = len(sigma)
percentages = [0.25, 0.50, 0.75]
k_values = [int(p * n_singular_values) for p in percentages]

# Create figure for original and reconstructed images
plt.figure(figsize=(15, 10))

# Plot original image
plt.subplot(2, 2, 1)
plt.imshow(A, cmap='gray')
plt.title('Original Image')
plt.axis('off')

# For each percentage, reconstruct and plot image
for i, (p, k) in enumerate(zip(percentages, k_values), 2):
    # Create truncated sigma
    sigma_k = np.copy(sigma)
    sigma_k[k:] = 0
    
    # Reconstruct image using k singular values
    A_k = U @ np.diag(sigma_k) @ Vt
    
    # Calculate Frobenius norm error
    error = np.linalg.norm(A - A_k, 'fro')
    
    # Plot reconstructed image
    plt.subplot(2, 2, i)
    plt.imshow(A_k, cmap='gray')
    plt.title(f'Reconstructed ({int(p*100)}% singular values)\nFrobenius Error: {error:.2f}')
    plt.axis('off')

plt.tight_layout()
plt.show()

# Print detailed error analysis
print("\nDetailed Error Analysis:")
print("-----------------------")
for p, k in zip(percentages, k_values):
    sigma_k = np.copy(sigma)
    sigma_k[k:] = 0
    A_k = U @ np.diag(sigma_k) @ Vt
    error = np.linalg.norm(A - A_k, 'fro')
    print(f"Using {k} singular values ({int(p*100)}%):")
    print(f"Frobenius norm error: {error:.2f}")
    print(f"Relative error: {error/np.linalg.norm(A, 'fro')*100:.2f}%\n")
