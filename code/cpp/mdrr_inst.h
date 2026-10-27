// MDRR 实例读入（与 prior/.../cpp/rost4.cpp 的 Inst::load 相同的文本格式；每站一行 "δ lead K"，K 可以不同）。
#pragma once
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef unsigned long long ull;

struct Type { int d, L, v0, flex; vector<int> elig; int q; vector<int> C; double post_h, back_h;
              vector<vector<pair<int,int>>> step; vector<int> pre; };
struct Site { double delta; int lead, K; };
struct ToolT { double ins, unin, maint; };

struct Inst {
    string name; int nk, ns, nt; vector<Type> T; vector<Site> S; vector<ToolT> TT;
    vector<int> N0, M0; vector<vector<int>> eidx; vector<int> nelig;
    // 读入：逐项与 rost4 相同
    void load(const string& path) {
        ifstream f(path); if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
        getline(f, name); f >> nk >> ns >> nt; T.resize(nk);
        for (int k = 0; k < nk; k++) {
            Type& t = T[k]; int ne; f >> t.d >> t.L >> t.v0 >> t.flex >> ne; t.elig.resize(ne);
            for (auto& e : t.elig) f >> e;
            int nc; f >> t.q >> nc; t.C.resize(nc); for (auto& c : t.C) f >> c;
            f >> t.post_h >> t.back_h;
            t.step.assign(t.L, vector<pair<int,int>>(ne));
            for (int i = 0; i < t.L; i++) for (int e = 0; e < ne; e++) f >> t.step[i][e].first >> t.step[i][e].second;
            t.pre.resize(ne); for (auto& p : t.pre) f >> p;
        }
        S.resize(ns); for (auto& s : S) f >> s.delta >> s.lead >> s.K;
        TT.resize(nt); for (auto& t : TT) f >> t.ins >> t.unin >> t.maint;
        N0.resize(nk); for (auto& x : N0) f >> x; M0.resize(nt); for (auto& x : M0) f >> x;
        if (!f) { fprintf(stderr, "bad instance file %s\n", path.c_str()); exit(1); }
        eidx.assign(nk, vector<int>(ns, -1)); nelig.assign(ns, 0);
        for (int k = 0; k < nk; k++) for (int e = 0; e < (int)T[k].elig.size(); e++) { eidx[k][T[k].elig[e]] = e; nelig[T[k].elig[e]]++; }
    }
};
