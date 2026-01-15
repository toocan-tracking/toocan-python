import numpy as np

def get_custom_kernel_3d():
    """
    Constructs a 3x3x3 boolean kernel with:
    - 8-connectivity in space (on the central time slice),
    - 2-connectivity in time (direct neighbors at t±1 in same spatial location).
    
    Returns:
        kernel (np.ndarray): 3D boolean array of shape (3, 3, 3)
    """
    kernel = np.zeros((3, 3, 3), dtype=bool)

    # 8-connected spatial neighbors in central time slice (t=0)
    kernel[1, :, :] = np.array([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]
    ])

    # Direct temporal neighbors at (x,y) position
    kernel[0, 1, 1] = 1  # t-1
    kernel[2, 1, 1] = 1  # t+1

    return kernel
