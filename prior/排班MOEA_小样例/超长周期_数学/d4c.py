import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[os.path.join(_h,'..','启发式验证'), os.path.join(_h,'..')]
import json, sys
from d4 import params, feasible
cases = [("L1-720",1,1,True,9),("L1-720",1,1,False,9),("L3-720",0,2,True,11),("L3-720",2,1,True,7),("L3-720",2,1,False,7),("L3-720",2,1,False,8),("L2-720",1,2,True,13),("L2-720",2,2,True,11)]
for name,r,K,early,N in cases:
    res = feasible(params(name,r), N, K, early, lim=6_000_000)
    print(json.dumps(dict(inst=name,type=r,K=K,提前维修=early,N=N,可行=res),ensure_ascii=False), flush=True)
