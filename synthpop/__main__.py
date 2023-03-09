import sys
import os

# include the synthpop main directory to path
pth = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, pth)

from synthpop.synthpop_main import main

main()
