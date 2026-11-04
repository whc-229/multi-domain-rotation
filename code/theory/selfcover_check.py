"""06 笔记 §3′：单机型自覆盖条件的穷举核对（全新启动 + 启动期提前大修，台数取速率下界 K·S/L）。
条件：对所有 1 ≤ i ≤ L−1、0 ≤ p ≤ L−2−i，有 g_p ≤ g_{i+p}。预期：满足条件的序列 K·S/L 台可行；[3,5,2,6] 不满足、不可行。
[4,6,5,4] 状态太多，未列入。"""
import sys, importlib.util
spec=importlib.util.spec_from_file_location("t",__file__.replace("selfcover_check.py","toy_connector.py"))
src=open(__file__.replace("selfcover_check.py","toy_connector.py")).read().split("for reset in")[0]
ns={}; exec(src,ns); feasible=ns["feasible"]
def cond(g):
    L=len(g); return all(g[p]<=g[i+p] for i in range(1,L) for p in range(0,L-1-i))
for g,K in [((4,4,5),3),((2,2,5),3),((3,5,2,6),4),((2,6,3),3)]:
    S=sum(g); L=len(g); N=K*S//L
    print(g,"K",K,"N=KS/L",N,"cond",cond(g),"fresh+reset feasible with N:",feasible(N,K,g,True), " no-reset:",feasible(N,K,g,False))
