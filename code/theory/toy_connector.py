"""06 笔记 §1 小例子：单机型 g = (2, 4)、两个岗位，全新启动时用状态图穷举判定 N 台能否永远排下去（不许 / 允许提前大修）。
预期输出：不许提前大修时 6 台不可行、7 台可行；允许提前大修时 6 台可行。"""
import itertools, sys
def feasible(N,K,g,reset):
    L=len(g)
    # robot state: (phase, remaining wait r>=0); ready if r==0. state = sorted tuple
    start=tuple(sorted([(0,0)]*N))
    def succ(st):
        st=list(st); ready=[i for i,(p,r) in enumerate(st) if r==0]
        out=set()
        for ch in itertools.combinations(ready,K):
            # each chosen robot: works now; options normal or reset (if p!=L-1)
            opts=[]
            for i in ch:
                p,_=st[i]; o=[((p+1)%L, g[p]-1)]
                if reset and p!=L-1: o.append((0,g[L-1]-1))
                opts.append(o)
            for combo in itertools.product(*opts):
                new=[]
                for i,(p,r) in enumerate(st):
                    if i in ch: continue
                    new.append((p,max(0,r-1)))
                new+= list(combo)
                out.add(tuple(sorted(new)))
        return out
    # find reachable states, then prune states without successors iteratively
    seen={start}; stack=[start]; E={}
    while stack:
        s=stack.pop(); E[s]=succ(s)
        for t in E[s]:
            if t not in seen: seen.add(t); stack.append(t)
    alive=set(seen)
    changed=True
    while changed:
        changed=False
        for s in list(alive):
            if not (E[s] & alive): alive.discard(s); changed=True
    return start in alive
for reset in (False,True):
    for N in (6,7,8):
        print("reset",reset,"N",N,feasible(N,2,(2,4),reset))
