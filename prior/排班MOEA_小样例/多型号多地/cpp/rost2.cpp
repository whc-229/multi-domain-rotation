// 多型号 × 多地 × 共享工具池的派工仿真、无限期证书（Brent 判环）、定台数、(μ+λ)-ES。
// 与 Python 版 certsim.py 逐行对应；实例由 export_inst.py 导出（所有取整后的步长表由 Python 预先算好）。
// 用法：
//   rost eval  <inst> <Tmax> <certify0/1> <Nk...> | <Mj...> | <w...>
//   rost es    <inst> <mode finite|cert> <H> <gens> <mu> <lam> <seed> <Tcert> <Tverify> [gp]
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef unsigned long long ull;

static const int CL = 40;   // 截断，与 Python 的 C 一致

struct Type { int d, L, v0, flex; vector<int> elig; int q; vector<int> C; double post_h, back_h;
              vector<vector<pair<int,int>>> step; vector<int> pre; };
struct Site { double delta; int lead, K; };
struct ToolT { double ins, unin, maint; };

struct Inst {
    string name; int nk, ns, nt; vector<Type> T; vector<Site> S; vector<ToolT> TT;
    vector<int> N0, M0; vector<vector<int>> eidx; vector<int> nelig;
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
        eidx.assign(nk, vector<int>(ns, -1)); nelig.assign(ns, 0);
        for (int k = 0; k < nk; k++) for (int e = 0; e < (int)T[k].elig.size(); e++) { eidx[k][T[k].elig[e]] = e; nelig[T[k].elig[e]]++; }
    }
};

// ---------------- 策略：线性打分（6 维机器人 + 2 维工具）或 GP 树（机器人打分）+ 线性工具打分 ----------------
struct GPNode { int op; double c; int a, b; };   // op: 0..5 特征, 6 常数, 7 +, 8 -, 9 *, 10 保护除, 11 max, 12 min
struct Policy {
    vector<double> w;            // 至少 8 维
    vector<GPNode> tree;         // 空则用线性
    double eval_tree(int n, const double* f) const {
        const GPNode& g = tree[n];
        if (g.op < 6) return f[g.op];
        if (g.op == 6) return g.c;
        double x = eval_tree(g.a, f), y = eval_tree(g.b, f);
        switch (g.op) {
            case 7: return x + y; case 8: return x - y; case 9: return x * y;
            case 10: return fabs(y) < 1e-9 ? 1.0 : x / y; case 11: return max(x, y); default: return min(x, y);
        }
    }
};

struct Res { ll gaps; int cert; ll period, t_end; vector<ll> why_r, why_t; ll enter = -1; };
static inline pair<ull,ull> hsh(const vector<ll>& v) { ull a = 1469598103934665603ULL, b = 0x9E3779B97F4A7C15ULL;
  for (ll x : v) { ull y = (ull)x; a ^= y + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2); a *= 1099511628211ULL; b = (b ^ (y * 0xff51afd7ed558ccdULL)) * 0xc4ceb9fe1a85ec53ULL; b ^= b >> 33; }
  return {a, b}; }
struct PH { size_t operator()(const pair<ull,ull>& p) const { return p.first ^ (p.second * 31); } };

static inline int postu(const Inst& I, int k, int s, int j, double unin) {
    const Type& t = I.T[k];
    return (int)ceil((t.post_h + t.back_h + I.S[s].delta + unin + I.TT[j].maint) / 2.0 - 1e-9);
}

Res run(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll Tmax, bool certify, ll stop_at = -1, vector<ll>* snap = nullptr) {
    int nk = I.nk, ns = I.ns, nt = I.nt;
    vector<int> rk, rph; vector<ll> rv;
    for (int k = 0; k < nk; k++) for (int c = 0; c < Nk[k]; c++) { rk.push_back(k); rph.push_back(0); rv.push_back(I.T[k].v0); }
    int nr = rk.size();
    bool use_tools = nt > 0;
    vector<vector<ll>> tf(nt);
    for (int j = 0; j < nt; j++) tf[j].assign(Mj[j], -1000000);
    int c0 = 0; for (int k = 0; k < nk; k++) c0 = max(c0, I.T[k].v0);
    vector<ll> fr(ns); for (int s = 0; s < ns; s++) fr[s] = c0 + I.S[s].lead;
    vector<vector<ll>> D(nk); vector<ll> tot(nk, 0);
    for (int k = 0; k < nk; k++) D[k].assign(I.T[k].L, 0);
    Res R; R.gaps = 0; R.cert = -1; R.period = -1; R.why_r.assign(nk, 0); R.why_t.assign(nt, 0);
    // Brent
    vector<ll> cur; unordered_map<pair<ull,ull>, ll, PH> seen;
    vector<int> due; vector<pair<double,int>> cand; vector<tuple<double,ll,int,int>> opts;
    ll t = *min_element(fr.begin(), fr.end());
    const bool lin = P.tree.empty();
    while (t < Tmax) {
        if (certify) {
            cur.clear();
            for (int r = 0; r < nr; r++) { cur.push_back(rk[r]); cur.push_back(rph[r]); cur.push_back(max(rv[r] - t, (ll)-CL)); }
            for (int s = 0; s < ns; s++) cur.push_back(fr[s] - t);
            for (int k = 0; k < nk; k++) for (int i = 0; i < I.T[k].L; i++) cur.push_back(tot[k] - (ll)I.T[k].L * D[k][i]);
            for (int j = 0; j < nt; j++) {
                size_t st0 = cur.size();
                for (ll f : tf[j]) cur.push_back(max(f - t, (ll)-CL));
                sort(cur.begin() + st0, cur.end());
            }
            if (stop_at >= 0) { if (t >= stop_at) { *snap = cur; R.t_end = t; return R; } }
            else {
                auto h = hsh(cur);
                auto it = seen.find(h);
                if (it != seen.end()) { R.cert = (R.gaps == 0) ? 1 : 0; R.period = t - it->second; R.enter = it->second; R.t_end = t; return R; }
                seen.emplace(h, t);
            }
        }
        due.clear();
        for (int s = 0; s < ns; s++) if (fr[s] <= t) due.push_back(s);
        stable_sort(due.begin(), due.end(), [&](int a, int b) { return I.nelig[a] < I.nelig[b]; });
        for (int s : due) {
            cand.clear();
            for (int r = 0; r < nr; r++) {
                int k = rk[r], e = I.eidx[k][s];
                if (e < 0 || rv[r] + I.S[s].lead > t) continue;
                int i = rph[r]; const Type& T = I.T[k];
                double f0 = (double)max(rv[r] - t, (ll)-CL);
                double nx = (double)T.step[i][e].second;
                double deficit = (double)tot[k] / (double)T.L - (double)D[k][i];
                double sc;
                if (lin) {
                    sc = P.w[0] * f0;
                    sc = sc + P.w[1] * nx / (double)T.d;
                    sc = sc + P.w[2] * (double)T.flex;
                    sc = sc - P.w[3] * deficit;
                    sc = sc + P.w[4] * (double)T.d;
                    sc = sc + P.w[5] * (double)i / (double)T.L;
                } else {
                    double fe[6] = {f0, nx / (double)T.d, (double)T.flex, deficit, (double)T.d, (double)i / (double)T.L};
                    sc = P.eval_tree(0, fe);
                    if (!isfinite(sc)) sc = 1e18;
                }
                cand.push_back({sc, r});
            }
            sort(cand.begin(), cand.end());
            bool done = false;
            for (auto& cr : cand) {
                int r = cr.second, k = rk[r], i = rph[r], e = I.eidx[k][s]; const Type& T = I.T[k];
                if (use_tools) {
                    ll need = t - T.pre[e];
                    opts.clear();
                    for (int j : T.C) for (int idx = 0; idx < (int)tf[j].size(); idx++) {
                        ll f = tf[j][idx];
                        if (f <= need) { ll rel = max(f - t, (ll)-CL); opts.emplace_back(P.w[6] * I.TT[j].maint + P.w[7] * (double)rel, rel, j, idx); }
                    }
                    if ((int)opts.size() < T.q) continue;
                    partial_sort(opts.begin(), opts.begin() + T.q, opts.end());
                    double unin = 0; for (int c = 0; c < T.q; c++) unin += I.TT[get<2>(opts[c])].unin;
                    for (int c = 0; c < T.q; c++) { int j = get<2>(opts[c]), idx = get<3>(opts[c]); tf[j][idx] = t + T.d + postu(I, k, s, j, unin); }
                }
                int ld = T.step[i][e].first, nx = T.step[i][e].second;
                rv[r] = t - ld + nx; rph[r] = (i + 1) % T.L;
                D[k][i]++; tot[k]++;
                fr[s] = t + T.d; done = true; break;
            }
            if (!done) {
                R.gaps++;
                if (!cand.empty()) { for (int j : I.T[rk[cand[0].second]].C) R.why_t[j]++; }
                else { for (int k = 0; k < nk; k++) if (I.eidx[k][s] >= 0) R.why_r[k]++; }
                fr[s] = t + 1;
            }
        }
        t = *min_element(fr.begin(), fr.end());
    }
    R.t_end = t; return R;
}

bool exact_check(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll enter, ll period) {
    vector<ll> a, b;
    Res r1 = run(I, Nk, Mj, P, LLONG_MAX / 4, true, enter, &a);
    Res r2 = run(I, Nk, Mj, P, LLONG_MAX / 4, true, enter + period, &b);
    return r1.t_end == enter && r2.t_end == enter + period && a == b;
}

struct SizeRes { int R; vector<int> Nk, Mj; string status; ll period; };

SizeRes size_res(const Inst& I, const Policy& P, ll H, const string& mode, ll Tcert, int max_add = 60) {
    vector<int> Nk = I.N0, Mj = I.M0; int nk = I.nk, nt = I.nt, nres = nk + nt;
    auto bump = [&](int idx, int d) { if (idx < nk) Nk[idx] += d; else Mj[idx - nk] += d; };
    auto tot = [&]() { int s = 0; for (int x : Nk) s += x; for (int x : Mj) s += x; return s; };
    Res r = run(I, Nk, Mj, P, H, false); int added = 0;
    while (r.gaps > 0 && added < max_add) {
        vector<ll> sig = r.why_r; sig.insert(sig.end(), r.why_t.begin(), r.why_t.end());
        ll bg = LLONG_MAX, bsig = 0; int bi = -1;
        for (int idx = 0; idx < nres; idx++) {
            bump(idx, 1); ll g = run(I, Nk, Mj, P, H, false).gaps; bump(idx, -1);
            if (bi < 0 || g < bg || (g == bg && sig[idx] > bsig)) { bg = g; bi = idx; bsig = sig[idx]; }
        }
        bump(bi, 1); added++; r = run(I, Nk, Mj, P, H, false);
    }
    if (r.gaps > 0) return {-1, Nk, Mj, "fail", -1};
    bool improved = true;
    while (improved) {
        improved = false;
        for (int idx = 0; idx < nres; idx++) {
            int& v = idx < nk ? Nk[idx] : Mj[idx - nk];
            if (v == 0) continue;
            v--;
            if (run(I, Nk, Mj, P, H, false).gaps == 0) { improved = true; break; }
            v++;
        }
    }
    if (mode == "finite") return {tot(), Nk, Mj, "finite", -1};
    for (int a = 0; a < max_add; a++) {
        Res c = run(I, Nk, Mj, P, Tcert, true);
        if (c.cert == 1) return {tot(), Nk, Mj, "cert", c.period};
        if (c.cert == -1 && c.gaps == 0) return {tot(), Nk, Mj, "uncert", -1};
        vector<ll> sig = c.why_r; sig.insert(sig.end(), c.why_t.begin(), c.why_t.end());
        int bi = max_element(sig.begin(), sig.end()) - sig.begin(); bump(bi, 1);
    }
    return {-1, Nk, Mj, "fail", -1};
}

// ---------------- GP 树随机生成与变异 ----------------
struct GP {
    mt19937_64& rng; GP(mt19937_64& r) : rng(r) {}
    int rand_int(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }
    double rand_u() { return uniform_real_distribution<double>(0, 1)(rng); }
    int grow(vector<GPNode>& t, int depth, int maxd) {
        int id = t.size(); t.push_back({0, 0, -1, -1});
        if (depth >= maxd || (depth > 1 && rand_u() < 0.3)) {
            if (rand_u() < 0.8) t[id] = {rand_int(0, 5), 0, -1, -1};
            else t[id] = {6, uniform_real_distribution<double>(-2, 2)(rng), -1, -1};
        } else {
            int op = rand_int(7, 12); int a = grow(t, depth + 1, maxd); int b = grow(t, depth + 1, maxd);
            t[id] = {op, 0, a, b};
        }
        return id;
    }
    // 子树抽取：把以 n 为根的子树复制到 out（重新编号），返回新根
    int copy_sub(const vector<GPNode>& t, int n, vector<GPNode>& out) {
        int id = out.size(); out.push_back(t[n]);
        if (t[n].op >= 7) { int a = copy_sub(t, t[n].a, out); int b = copy_sub(t, t[n].b, out); out[id].a = a; out[id].b = b; }
        return id;
    }
    int depth(const vector<GPNode>& t, int n) { return t[n].op >= 7 ? 1 + max(depth(t, t[n].a), depth(t, t[n].b)) : 1; }
    // 把 t 中节点 target 替换为 sub 子树，返回新树
    vector<GPNode> replace(const vector<GPNode>& t, int target, const vector<GPNode>& sub) {
        vector<GPNode> out; function<int(int)> rec = [&](int n) -> int {
            if (n == target) return copy_sub(sub, 0, out);
            int id = out.size(); out.push_back(t[n]);
            if (t[n].op >= 7) { int a = rec(t[n].a); int b = rec(t[n].b); out[id].a = a; out[id].b = b; }
            return id;
        };
        rec(0); return out;
    }
    vector<GPNode> mutate(const vector<GPNode>& t) {
        int target = rand_int(0, (int)t.size() - 1);
        vector<GPNode> sub; grow(sub, 0, 3);
        auto out = replace(t, target, sub);
        if (depth(out, 0) > 7) return t;
        return out;
    }
    vector<GPNode> crossover(const vector<GPNode>& a, const vector<GPNode>& b) {
        int ta = rand_int(0, (int)a.size() - 1), tb = rand_int(0, (int)b.size() - 1);
        vector<GPNode> sub; copy_sub(b, tb, sub);
        auto out = replace(a, ta, sub);
        if (depth(out, 0) > 7) return a;
        return out;
    }
    string str(const vector<GPNode>& t, int n) {
        static const char* nm[6] = {"wait", "nextocc", "flex", "deficit", "duty", "phase"};
        const GPNode& g = t[n];
        if (g.op < 6) return nm[g.op];
        if (g.op == 6) { char b[32]; snprintf(b, 32, "%.3g", g.c); return b; }
        static const char* on[6] = {"+", "-", "*", "/", "max", "min"};
        if (g.op >= 11) return string(on[g.op - 7]) + "(" + str(t, g.a) + "," + str(t, g.b) + ")";
        return "(" + str(t, g.a) + on[g.op - 7] + str(t, g.b) + ")";
    }
};

string vec_s(const vector<int>& v) { string s = "["; for (size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + to_string(v[i]); return s + "]"; }
string dvec_s(const vector<double>& v) { string s = "["; char b[32]; for (size_t i = 0; i < v.size(); i++) { snprintf(b, 32, "%s%.4f", i ? "," : "", v[i]); s += b; } return s + "]"; }

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    string cmd = argv[1]; Inst I; I.load(argv[2]);
    if (cmd == "eval") {
        ll Tmax = atoll(argv[3]); bool cert = atoi(argv[4]);
        vector<int> Nk, Mj; vector<double> w; int part = 0;
        for (int a = 5; a < argc; a++) {
            string s = argv[a]; if (s == "|") { part++; continue; }
            if (part == 0) Nk.push_back(atoi(s.c_str())); else if (part == 1) Mj.push_back(atoi(s.c_str())); else w.push_back(atof(s.c_str()));
        }
        while (w.size() < 8) w.push_back(0);
        Policy P; P.w = w;
        Res r = run(I, Nk, Mj, P, Tmax, cert);
        int exact = (r.cert == 1) ? (exact_check(I, Nk, Mj, P, r.enter, r.period) ? 1 : 0) : -1;
        printf("{\"gaps\": %lld, \"cert\": %d, \"period\": %lld, \"t_end\": %lld, \"exact\": %d}\n", r.gaps, r.cert, r.period, r.t_end, exact);
        return 0;
    }
    if (cmd == "es") {
        string mode = argv[3]; ll H = atoll(argv[4]); int gens = atoi(argv[5]), mu = atoi(argv[6]), lam = atoi(argv[7]);
        ull seed = strtoull(argv[8], 0, 10); ll Tcert = atoll(argv[9]), Tver = atoll(argv[10]);
        bool gp = argc > 11 && string(argv[11]) == "gp";
        mt19937_64 rng(seed); normal_distribution<double> N01(0, 1); GP G(rng);
        auto t0 = chrono::steady_clock::now();
        struct Ind { Policy P; double f; SizeRes s; };
        auto fit = [&](Policy& P) { SizeRes s = size_res(I, P, H, mode, Tcert); double f = s.R < 0 ? 999 : s.R + (s.status == "uncert" ? 0.5 : 0); return make_pair(f, s); };
        vector<Ind> pop;
        vector<vector<double>> base = {{1, 0, 0, 0, 0, 0, 1, 0}, {0.1, 0, 10, 1, 0, 0, 1, 0.1}, {1, 0.5, 1, 0.5, 0, 0, 1, -0.1}};
        for (int a = 0; a < mu + lam; a++) {
            Policy P; P.w.assign(8, 0);
            if (a < 3) P.w = base[a]; else for (auto& x : P.w) x = N01(rng);
            if (gp) {
                if (a == 0) P.tree = {{0, 0, -1, -1}};                                                // FIFO: wait
                else if (a == 1) P.tree = {{8, 0, 1, 2}, {9, 0, 3, 4}, {6, 10, -1, -1}, {2, 0, -1, -1}, {3, 0, -1, -1}}; // 10*flex - deficit 近似 PB
                else { G.grow(P.tree, 0, 4); }
                P.w[6] = 1; P.w[7] = (a == 1) ? 0.1 : 0;
            }
            auto fs = fit(P); pop.push_back({P, fs.first, fs.second});
        }
        sort(pop.begin(), pop.end(), [](const Ind& a, const Ind& b) { return a.f < b.f; }); pop.resize(mu);
        vector<double> hist = {pop[0].f};
        for (int g = 0; g < gens; g++) {
            vector<Ind> kids;
            for (int c = 0; c < lam; c++) {
                const Ind& pa = pop[uniform_int_distribution<int>(0, mu - 1)(rng)];
                Policy P = pa.P;
                if (!gp) { for (auto& x : P.w) x = x + 0.5 * N01(rng) * (fabs(x) + 0.2); }
                else {
                    const Ind& pb = pop[uniform_int_distribution<int>(0, mu - 1)(rng)];
                    if (uniform_real_distribution<double>(0, 1)(rng) < 0.5) P.tree = G.crossover(pa.P.tree, pb.P.tree);
                    else P.tree = G.mutate(pa.P.tree);
                    for (int q = 6; q < 8; q++) P.w[q] = P.w[q] + 0.3 * N01(rng) * (fabs(P.w[q]) + 0.2);
                }
                auto fs = fit(P); kids.push_back({P, fs.first, fs.second});
            }
            for (auto& k : kids) pop.push_back(k);
            stable_sort(pop.begin(), pop.end(), [](const Ind& a, const Ind& b) { return a.f < b.f; }); pop.resize(mu);
            hist.push_back(pop[0].f);
        }
        double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        const Ind& b = pop[0];
        Res v = (b.s.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, b.s.Nk, b.s.Mj, b.P, Tver, true);
        int exact = -1;
        if (v.cert == 1) exact = exact_check(I, b.s.Nk, b.s.Mj, b.P, v.enter, v.period) ? 1 : 0;
        string hs = "["; char buf[32]; for (size_t i = 0; i < hist.size(); i++) { snprintf(buf, 32, "%s%.1f", i ? "," : "", hist[i]); hs += buf; } hs += "]";
        printf("{\"inst\": \"%s\", \"mode\": \"%s\", \"H\": %lld, \"gp\": %d, \"seed\": %llu, \"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, "
               "\"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"sec\": %.1f, \"hist\": %s, \"w\": %s, \"tree\": \"%s\"}\n",
               I.name.c_str(), mode.c_str(), H, gp ? 1 : 0, seed, b.s.R, b.s.status.c_str(), vec_s(b.s.Nk).c_str(), vec_s(b.s.Mj).c_str(),
               v.gaps, v.cert, v.period, exact, sec, hs.c_str(), dvec_s(b.P.w).c_str(), gp ? G.str(b.P.tree, 0).c_str() : "");
        return 0;
    }
    return 1;
}
