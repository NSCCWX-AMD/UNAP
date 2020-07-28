# -*- coding: utf-8 -*-
"""
SCons build utilities
----------------------------

Collection of miscellaneous utilities that don't fit elsewhere.
"""

import os
import platform


def ostype():
    """Return the operating system type"""

    if os.name == 'nt':
        return 'windows'
    else:
        return os.uname()[0].lower()
