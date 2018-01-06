# updates system path with the project's modules if imported in a notebook

import sys
import pathlib
sys.path.insert(0,str(pathlib.Path(__file__).parents[1]))
