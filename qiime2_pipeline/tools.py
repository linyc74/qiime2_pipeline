import os
from typing import List


def get_files(
        source: str = '.',
        startswith: str = '',
        endswith: str = '',
        isfullpath: bool = False) -> List[str]:
    """
    Get all files that start with or end with some strings in the source folder

    Args:
        startswith

        endswith

        source: path-like
            The source directory

        isfullpath
            If True, return the full file paths in the list

    Returns:
        A list of file names, or file paths
        If no files found, return an empty list []
    """
    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [f for f in files if (f.startswith(s) and f.endswith(e))]

    if isfullpath:
        ret = [os.path.join(source, f) for f in ret]

    if ret:
        # Sort the file list so that the order will be
        #     consistent across different OS platforms
        ret.sort()
    return ret


def rev_comp(seq: str) -> str:
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N',
        'M': 'K',  # M = A C
        'K': 'M',  # K = G T
        'R': 'Y',  # R = A G
        'Y': 'R',  # Y = C T
        'S': 'S',  # S = C G
        'W': 'W',  # W = A T
        'B': 'V',  # B = C G T
        'V': 'B',  # V = A C G
        'D': 'H',  # D = A G T
        'H': 'D',  # H = A C T

        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'n': 'n',
        'm': 'k',
        'k': 'm',
        'r': 'y',
        'y': 'r',
        's': 's',
        'w': 'w',
        'b': 'v',
        'v': 'b',
        'd': 'h',
        'h': 'd'
    }
    return ''.join([comp[base] for base in seq[::-1]])
