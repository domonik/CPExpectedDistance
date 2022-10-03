from CPExpectedDistance.c_expected_distance import cp_expected_distance
import RNA


def expected_distance(fc: RNA.fold_compound):
    """Python wrapper for Clote-Ponty expected distance calculation

    Args:
        fc (RNA.fold_compound): The ViennaRNA Fold compound

    Returns (np.ndarray): Numpy Array containing Expected Distances

    """
    # fc.this is the pointer to the Swig interface object from which we have access to the underlying
    # ViennaRNA C fold compound
    exp_d = cp_expected_distance(fc.this)
    return exp_d
