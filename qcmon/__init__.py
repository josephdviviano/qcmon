"""
gotta qc 'em all.
"""
import os, glob

modules = glob.glob(os.path.dirname(__file__)+"/*.py")
__all__ = [ os.path.basename(f)[:-3] for f in modules]

from . import *
