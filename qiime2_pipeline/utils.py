import os
from typing import Optional


def edit_fpath(
        fpath: str,
        old_suffix: str = '',
        new_suffix: str = '',
        dstdir: Optional[str] = None) -> str:

    f = os.path.basename(fpath)
    if old_suffix != '':
        f = f[:-len(old_suffix)]  # strip old suffix
    f += new_suffix

    if dstdir is None:
        dstdir = os.path.dirname(fpath)

    return f'{dstdir}/{f}'


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:

    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1
