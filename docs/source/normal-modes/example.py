import numpy as np

def example(arr, axis, indices, keepdims=False):
    """
    :param arr: _description_
    :type arr: _type_
    :param axis: _description_
    :type axis: _type_
    :param indices: _description_
    :type indices: _type_
    :param keepdims: _description_, defaults to False
    :type keepdims: bool, optional
    :return: _description_
    :rtype: _type_
    """
    result = np.take(arr, indices=indices, axis=axis)
    
    # If indices is a scalar, add back the singleton dimension along the axis
    if np.isscalar(indices) and keepdims:
        result = np.expand_dims(result, axis=axis)
    
    return result
