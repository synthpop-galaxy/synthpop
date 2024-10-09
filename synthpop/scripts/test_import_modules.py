import sys
import os

pth = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, pth)
import modules

for i in os.listdir("../../modules"):
    if not os.path.isdir(f"../modules/{i}"): continue
    if i.startswith("__"): continue
    for j in os.listdir(f"../modules/{i}"):
        if j.startswith("__"): continue
        if not j.endswith(".py"): continue
        try:
            exec(f"from modules.{i} import {j[:-3]}")
        except (ValueError, ImportError) as e:
            print(f"can not import {j} from {i}")
            print(e)
