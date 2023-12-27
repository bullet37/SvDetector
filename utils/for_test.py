def get_median(elements):
    elements = sorted(elements)
    middle = len(elements) // 2
    return elements[middle]


def get_percentile(vec, p):
    idx = int(len(vec) * p)
    return sorted(vec)[idx]