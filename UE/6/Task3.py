import numpy as np
import matplotlib.pyplot as plt

def trig_interpolation(y):
    n = len(y)
    x = 2 * np.pi * np.arange(n) / n  # Interpolation points
    
    # Number of coefficients (n must be odd for this implementation)
    # For even n, we would have an asymmetric number of positive and negative frequencies,
    # which complicates the coefficient extraction. With odd n, we have m positive and m
    # negative frequencies plus the zero frequency, making the math simpler
    m = (n - 1) // 2
    
    # Compute coefficients using FFT
    # Divide by n to normalize the FFT coefficients according to the standard Fourier series formula
    # Without this normalization, the coefficients would be n times too large
    # This matches the continuous Fourier transform definition: 1/T ∫f(t)e^(-2πikt/T)dt
    c = np.fft.fft(y) / n
    
    # Extract real coefficients ak and bk
    ak = np.zeros(m + 1)
    bk = np.zeros(m + 1)
    
    # a0 is real part of c[0]
    ak[0] = np.real(c[0])
    
    # The factor of 2 appears because the FFT combines positive and negative frequencies.
    # When we convert from complex Fourier coefficients (c_k) to real Fourier coefficients (ak, bk),
    # we need to account for both c_k and c_{-k} which are complex conjugates.
    # This leads to: ak = c_k + c_{-k} = 2*Re(c_k) and bk = i(c_k - c_{-k}) = -2*Im(c_k)
    for k in range(1, m + 1):
        ak[k] = 2 * np.real(c[k])
        bk[k] = -2 * np.imag(c[k])
    
    # Create fine grid for plotting
    x_fine = np.linspace(0, 2*np.pi, 1000)
    y_interp = ak[0] / 2  # The a0 coefficient is divided by 2 in the Fourier series formula to avoid counting it twice when summing
    
    # Compute interpolant
    for k in range(1, m + 1):
        y_interp += ak[k] * np.cos(k * x_fine) + bk[k] * np.sin(k * x_fine)
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(x_fine, y_interp, 'b-', label='Interpolant')
    plt.plot(x, y, 'ro', label='Data points')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Trigonometric Interpolation')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return ak, bk

# Example usage
if __name__ == "__main__":
    # Test with sample data (e.g., sampling from sin(x))
    n = 11  # odd number of points
    x_sample = 2 * np.pi * np.arange(n) / n
    y_sample = 2*np.sin(2*x_sample) + 0.5*np.cos(3*x_sample) + 0.3*np.sin(4*x_sample) - 0.7*np.cos(x_sample) + 1.5
    
    ak, bk = trig_interpolation(y_sample)
    print("ak coefficients:", ak)
    print("bk coefficients:", bk)
