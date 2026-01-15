import numpy as np
from toocan.detection_spreading.detection import run_labeling_with_c

def test_run_labeling_with_c():
    print("\n?? Running test: `test_run_labeling_with_c`")

    shape = (5, 10, 10)
    volume = np.zeros(shape, dtype=bool)

    # Blob 1 at z=1
    volume[1, 1:4, 1:4] = True

    # Blob 2 at z=3 (spatially and temporally far)
    volume[3, 6:9, 6:9] = True

    coords = np.argwhere(volume)
    labels = run_labeling_with_c(coords, shape)
    n_objs = len(np.unique(labels))

    print("  ? Unique labels with two blobs:", np.unique(labels))
    assert n_objs == 2, f"Expected 2 objects, got {n_objs}"
    print("? Test passed: 2 blobs correctly separated.")


if __name__ == "__main__":
    test_run_labeling_with_c()