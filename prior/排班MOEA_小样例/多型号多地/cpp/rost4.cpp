// rost4：在 rost3 基础上加"提前大修"决策（reset）和"按当前相位人数"的相位特征（curdef）。
// 两个开关都关闭时与 rost3 完全相同。
// reset：机器人在相位 i ≠ L−1 被派去值岗时，策略决定值完这一班是否直接做最高级维修（占用按相位 L−1 计，回到相位 0）。
//   打分 = w8 + w9·(本相位人数占比 − 1/L) + w10·(相位 0 人数占比 − 1/L) + w11·i/(L−1) + w12·min(t/Tw, 1)，> 0 则提前大修。
//   用到 t 的特征只在 Tw 之前变化，所以判环从 t ≥ Tw 才开始记状态。hard：只在 t < Tw 时允许提前大修（只在启动期错峰）。
// curdef：派工打分里的"相位缺额"改为 当前在相位 i 的同型台数 − N_k/L（有界，提前大修时累计缺额会发散，状态不再有限）。
// stagger：错峰首修（由精确解读出的结构）：每型第 c 台机器人分到第 c mod L 组，第 g 组在首个超周期里第 g+1 班后就做最高级维修，之后照常。
// 用法：
//   rost4 eval <inst> <Tmax> <certify0/1> <Nk...> | <Mj...> | <w...> [| opts]
//   rost4 es   <inst> <mode finite|cert|certs> <H> <gens> <mu> <lam> <seed> <Tcert> <Tverify> [gp] [reset] [curdef] [Tw=格数]
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
    vector<double> w;            // 至少 8 维；reset 时 13 维
    vector<GPNode> tree;         // 空则用线性
    bool reset = false, curdef = false, hard = false, stagger = false; ll Tw = 0;
    vector<vector<int>> quota; vector<vector<ll>> qfrom;
    vector<vector<vector<pair<ll,double>>>> planf; vector<double> planN;
    vector<vector<double>> lpz;   // lpz[k][s]：LP 解里 k 型承担 s 站值岗的比例（LP 引导派工，权重 w13）   // 按比例的计划：LP 流量与 LP 台数，运行时按实际台数换算
    vector<vector<vector<pair<ll,int>>>> plan;   // plan[k][i] = [(t0, c)...] 按 t0 排序：t ≥ t0 后还有余量就提前大修   // quota[k][i]：k 型在相位 i 值岗后提前大修的总次数（来自全新启动 LP/ILP），只在启动期用
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

struct Res { ll gaps; int cert; ll period, t_end; vector<ll> why_r, why_t; ll enter = -1;
             ll resets = 0, duties = 0, cyc_resets = -1, cyc_duties = -1; };
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
    vector<int> grp; vector<char> firstc; vector<vector<int>> qleft = P.quota; auto pl = P.plan;
    if (!P.planf.empty()) {   // 按实际台数换算：累计流量 × Nk/N_LP 四舍五入
        pl.assign(nk, {});
        for (int k = 0; k < nk; k++) {
            pl[k].assign(I.T[k].L, {});
            double sc = (P.planN[k] > 1e-9) ? (double)Nk[k] / P.planN[k] : 0.0;
            for (int i = 0; i < I.T[k].L; i++) { double F = 0; int prev = 0;
                for (auto& e2 : P.planf[k][i]) { F += e2.second * sc; int c = (int)llround(F) - prev; if (c > 0) { pl[k][i].push_back({e2.first, c}); prev += c; } } }
        }
    }
    for (int k = 0; k < nk; k++) for (int c = 0; c < Nk[k]; c++) { rk.push_back(k); rph.push_back(0); rv.push_back(I.T[k].v0); grp.push_back(c % I.T[k].L); firstc.push_back(1); }
    int nr = rk.size();
    bool use_tools = nt > 0;
    vector<vector<ll>> tf(nt);
    for (int j = 0; j < nt; j++) tf[j].assign(Mj[j], -1000000);
    int c0 = 0; for (int k = 0; k < nk; k++) c0 = max(c0, I.T[k].v0);
    vector<ll> fr(ns); for (int s = 0; s < ns; s++) fr[s] = c0 + I.S[s].lead;
    vector<vector<ll>> D(nk); vector<ll> tot(nk, 0);
    for (int k = 0; k < nk; k++) D[k].assign(I.T[k].L, 0);
    Res R; R.gaps = 0; R.cert = -1; R.period = -1; R.why_r.assign(nk, 0); R.why_t.assign(nt, 0);
    vector<vector<int>> cnt(nk); for (int k = 0; k < nk; k++) { cnt[k].assign(I.T[k].L, 0); cnt[k][0] = Nk[k]; }
    const bool use_tw = P.reset && P.Tw > 0;
    ll qmax = 0; for (auto& v : P.qfrom) for (ll x : v) qmax = max(qmax, x);
    for (auto& pk : P.plan) for (auto& pi : pk) for (auto& e2 : pi) qmax = max(qmax, e2.first);
    for (auto& pk : P.planf) for (auto& pi : pk) for (auto& e2 : pi) qmax = max(qmax, e2.first);   // 配额起始时刻之前策略依赖 t，不判环
    // 判环：状态第一次重复
    vector<ll> cur; unordered_map<pair<ull,ull>, ll, PH> seen; vector<ll> st_res, st_dut;
    vector<int> due; vector<pair<double,int>> cand; vector<tuple<double,ll,int,int>> opts;
    ll t = *min_element(fr.begin(), fr.end());
    const bool lin = P.tree.empty();
    while (t < Tmax) {
        if (certify) {
            cur.clear();
            for (int r = 0; r < nr; r++) { cur.push_back(rk[r]); cur.push_back(rph[r]); cur.push_back(max(rv[r] - t, (ll)-CL)); if (P.stagger) cur.push_back(firstc[r]); }
            for (auto& qk : qleft) for (int q : qk) cur.push_back(q);
            for (auto& pk : pl) for (auto& pi : pk) for (auto& e2 : pi) cur.push_back(e2.second);
            for (int s = 0; s < ns; s++) cur.push_back(fr[s] - t);
            if (!P.curdef) for (int k = 0; k < nk; k++) for (int i = 0; i < I.T[k].L; i++) cur.push_back(tot[k] - (ll)I.T[k].L * D[k][i]);
            for (int j = 0; j < nt; j++) {
                size_t st0 = cur.size();
                for (ll f : tf[j]) cur.push_back(max(f - t, (ll)-CL));
                sort(cur.begin() + st0, cur.end());
            }
            if (stop_at >= 0) { if (t >= stop_at) { *snap = cur; R.t_end = t; return R; } }
            else if ((!use_tw || t >= P.Tw) && t >= qmax) {
                auto h = hsh(cur);
                auto it = seen.find(h);
                if (it != seen.end()) {
                    ll idx = it->second, t1 = idx >> 20;   // 高位存时刻、低位存序号
                    ll q = idx & ((1LL << 20) - 1);
                    R.cert = (R.gaps == 0) ? 1 : 0; R.period = t - t1; R.enter = t1; R.t_end = t;
                    R.cyc_resets = R.resets - st_res[q]; R.cyc_duties = R.duties - st_dut[q]; return R;
                }
                if (st_res.size() >= (1u << 20)) { R.t_end = t; return R; }   // 状态表满：按"未判出"处理
                seen.emplace(h, (t << 20) | (ll)st_res.size()); st_res.push_back(R.resets); st_dut.push_back(R.duties);
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
                double deficit = P.curdef ? ((double)cnt[k][i] - (double)Nk[k] / (double)T.L)
                                          : ((double)tot[k] / (double)T.L - (double)D[k][i]);
                double sc;
                if (lin) {
                    sc = P.w[0] * f0;
                    sc = sc + P.w[1] * nx / (double)T.d;
                    sc = sc + P.w[2] * (double)T.flex;
                    sc = sc - P.w[3] * deficit;
                    sc = sc + P.w[4] * (double)T.d;
                    sc = sc + P.w[5] * (double)i / (double)T.L;
                    if (!P.lpz.empty()) sc = sc - P.w[13] * P.lpz[k][s];
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
                int ld = T.step[i][e].first, nx = T.step[i][e].second, nph = (i + 1) % T.L;
                bool planned = false;
                if (!pl.empty() && i != T.L - 1) {
                    for (auto& e2 : pl[k][i]) { if (e2.first > t) break; if (e2.second > 0) { e2.second--; planned = true; break; } }
                }
                if (planned) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; }
                else if (!pl.empty()) {}
                else if (!qleft.empty()) {   // 配额：k 型相位 i 还有配额就提前大修
                    if (i != T.L - 1 && qleft[k][i] > 0 && t >= P.qfrom[k][i]) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; qleft[k][i]--; }
                } else if (P.stagger) {   // 错峰首修：第 g 组在首个超周期的第 g+1 班后提前大修，之后照常
                    if (firstc[r]) {
                        if (i == grp[r] && i != T.L - 1) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; firstc[r] = 0; }
                        else if (i == T.L - 1) firstc[r] = 0;
                    }
                } else if (P.reset && i != T.L - 1 && Nk[k] > 0 && (!P.hard || t < P.Tw)) {
                    double inv = 1.0 / (double)T.L, nkd = (double)Nk[k];
                    double sc = P.w[8] + P.w[9] * ((double)cnt[k][i] / nkd - inv) + P.w[10] * ((double)cnt[k][0] / nkd - inv)
                              + P.w[11] * (double)i / (double)(T.L - 1) + P.w[12] * (P.Tw > 0 ? min((double)t / (double)P.Tw, 1.0) : 1.0);
                    if (sc > 0) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; }
                }
                rv[r] = t - ld + nx; rph[r] = nph;
                cnt[k][i]--; cnt[k][nph]++; R.duties++;
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

vector<int> g_N0, g_M0;   // LP 引导时的定台数起点
SizeRes size_res(const Inst& I, const Policy& P, ll H, const string& mode, ll Tcert, int max_add = 60) {
    vector<int> Nk = g_N0.empty() ? I.N0 : g_N0, Mj = g_M0.empty() ? I.M0 : g_M0; int nk = I.nk, nt = I.nt, nres = nk + nt;
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
    bool strict = (mode == "certs");
    for (int a = 0; a < max_add; a++) {
        Res c = run(I, Nk, Mj, P, Tcert, true);
        if (strict && c.cert == -1 && c.gaps == 0) c = run(I, Nk, Mj, P, 4 * Tcert, true);   // 判不出循环：加长到 4 倍再判
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

vector<vector<ll>> g_qfrom; vector<double> g_planN;
vector<vector<double>> parse_lpz(const Inst& I, const string& fn, vector<int>& N0, vector<int>& M0) {   // 行 "z k s 比例"；"#N0 …"；"#M0 …"
    vector<vector<double>> z(I.nk, vector<double>(I.ns, 0)); ifstream f(fn); string line;
    while (getline(f, line)) {
        if (line.rfind("#N0", 0) == 0) { stringstream ss(line.substr(3)); N0.assign(I.nk, 0); for (auto& x : N0) ss >> x; }
        else if (line.rfind("#M0", 0) == 0) { stringstream ss(line.substr(3)); M0.assign(I.nt, 0); for (auto& x : M0) ss >> x; }
        else { int k, s2; double v; if (sscanf(line.c_str(), "z %d %d %lf", &k, &s2, &v) == 3 && k < I.nk && s2 < I.ns) z[k][s2] = v; }
    }
    return z;
}
vector<vector<vector<pair<ll,double>>>> parse_planf(const Inst& I, const string& fn) {   // 首行 "#N n0 n1 …"，其余每行：k i t0 流量
    vector<vector<vector<pair<ll,double>>>> p(I.nk); for (int k = 0; k < I.nk; k++) p[k].assign(I.T[k].L, {});
    ifstream f(fn); string line; g_planN.assign(I.nk, 0);
    while (getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') { stringstream ss(line.substr(2)); for (int k = 0; k < I.nk; k++) ss >> g_planN[k]; continue; }
        int k, i; long long t0; double x; if (sscanf(line.c_str(), "%d %d %lld %lf", &k, &i, &t0, &x) == 4 && k < I.nk && i < I.T[k].L) p[k][i].push_back({t0, x});
    }
    for (auto& pk : p) for (auto& pi : pk) sort(pi.begin(), pi.end());
    return p;
}
vector<vector<vector<pair<ll,int>>>> parse_plan(const Inst& I, const string& fn) {   // 文件每行：k i t0 c
    vector<vector<vector<pair<ll,int>>>> p(I.nk); for (int k = 0; k < I.nk; k++) p[k].assign(I.T[k].L, {});
    ifstream f(fn); int k, i, c; long long t0;
    while (f >> k >> i >> t0 >> c) if (k < I.nk && i < I.T[k].L && c > 0) p[k][i].push_back({t0, c});
    for (auto& pk : p) for (auto& pi : pk) sort(pi.begin(), pi.end());
    return p;
}
vector<vector<int>> parse_quota(const Inst& I, const string& str) {   // "k:i:c" 或 "k:i:c:t0"（t0 起才用）
    vector<vector<int>> q(I.nk); g_qfrom.assign(I.nk, {}); for (int k = 0; k < I.nk; k++) { q[k].assign(I.T[k].L, 0); g_qfrom[k].assign(I.T[k].L, 0); }
    stringstream ss(str); string item;
    while (getline(ss, item, ';')) { int k, i, c; long long t0 = 0; int n = sscanf(item.c_str(), "%d:%d:%d:%lld", &k, &i, &c, &t0); if (n >= 3) { q[k][i] += c; g_qfrom[k][i] = t0; } }
    return q;
}
string vec_s(const vector<int>& v) { string s = "["; for (size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + to_string(v[i]); return s + "]"; }
string dvec_s(const vector<double>& v) { string s = "["; char b[32]; for (size_t i = 0; i < v.size(); i++) { snprintf(b, 32, "%s%.4f", i ? "," : "", v[i]); s += b; } return s + "]"; }

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    string cmd = argv[1]; Inst I; I.load(argv[2]);
    if (cmd == "eval") {
        ll Tmax = atoll(argv[3]); bool cert = atoi(argv[4]);
        vector<int> Nk, Mj; vector<double> w; int part = 0; Policy P;
        for (int a = 5; a < argc; a++) {
            string s = argv[a]; if (s == "|") { part++; continue; }
            if (part == 0) Nk.push_back(atoi(s.c_str())); else if (part == 1) Mj.push_back(atoi(s.c_str()));
            else if (part == 2) w.push_back(atof(s.c_str()));
            else { if (s == "reset") P.reset = true; else if (s == "curdef") P.curdef = true; else if (s == "hard") P.hard = true; else if (s == "stagger") { P.stagger = true; P.curdef = true; } else if (s.rfind("quota=", 0) == 0) { P.curdef = true; P.quota = parse_quota(I, s.substr(6)); P.qfrom = g_qfrom; } else if (s.rfind("plan=", 0) == 0) { P.curdef = true; P.plan = parse_plan(I, s.substr(5)); } else if (s.rfind("planf=", 0) == 0) { P.curdef = true; P.planf = parse_planf(I, s.substr(6)); P.planN = g_planN; } else if (s.rfind("Tw=", 0) == 0) P.Tw = atoll(s.c_str() + 3); }
        }
        if (P.reset) P.curdef = true;
        while (w.size() < 13) w.push_back(0);
        P.w = w;
        Res r = run(I, Nk, Mj, P, Tmax, cert);
        int exact = (r.cert == 1) ? (exact_check(I, Nk, Mj, P, r.enter, r.period) ? 1 : 0) : -1;
        printf("{\"gaps\": %lld, \"cert\": %d, \"period\": %lld, \"t_end\": %lld, \"exact\": %d, \"resets\": %lld, \"duties\": %lld, \"cyc_resets\": %lld, \"cyc_duties\": %lld}\n",
               r.gaps, r.cert, r.period, r.t_end, exact, r.resets, r.duties, r.cyc_resets, r.cyc_duties);
        return 0;
    }
    if (cmd == "es") {
        string mode = argv[3]; ll H = atoll(argv[4]); int gens = atoi(argv[5]), mu = atoi(argv[6]), lam = atoi(argv[7]);
        ull seed = strtoull(argv[8], 0, 10); ll Tcert = atoll(argv[9]), Tver = atoll(argv[10]);
        bool gp = false, o_reset = false, o_curdef = false, o_hard = false, o_pad = false, o_stag = false; ll o_Tw = 0; vector<vector<int>> o_quota; vector<vector<vector<pair<ll,int>>>> o_plan; vector<vector<vector<pair<ll,double>>>> o_planf; vector<vector<double>> o_lpz;
        for (int a = 11; a < argc; a++) { string s = argv[a];
            if (s == "gp") gp = true; else if (s == "reset") o_reset = true; else if (s == "curdef") o_curdef = true; else if (s == "hard") o_hard = true; else if (s == "pad") o_pad = true; else if (s == "stagger") o_stag = true; else if (s.rfind("quota=", 0) == 0) o_quota = parse_quota(I, s.substr(6)); else if (s.rfind("plan=", 0) == 0) o_plan = parse_plan(I, s.substr(5)); else if (s.rfind("planf=", 0) == 0) o_planf = parse_planf(I, s.substr(6)); else if (s.rfind("lpg=", 0) == 0) o_lpz = parse_lpz(I, s.substr(4), g_N0, g_M0);
            else if (s.rfind("Tw=", 0) == 0) o_Tw = atoll(s.c_str() + 3); }
        if (o_reset || o_stag || !o_quota.empty() || !o_plan.empty() || !o_planf.empty()) o_curdef = true;
        const int NW = !o_lpz.empty() ? 14 : ((o_reset || o_pad) ? 13 : 8);   // pad：维度和随机数流与 reset 相同，但从不提前大修（对照组）
        mt19937_64 rng(seed); normal_distribution<double> N01(0, 1); GP G(rng);
        auto t0 = chrono::steady_clock::now();
        struct Ind { Policy P; double f; SizeRes s; };
        auto fit = [&](Policy& P) { SizeRes s = size_res(I, P, H, mode, Tcert); double pen = (mode == "certs") ? 3.0 : 0.5;
            double f = s.R < 0 ? 999 : s.R + (s.status == "uncert" ? pen : 0); return make_pair(f, s); };
        vector<Ind> pop;
        vector<vector<double>> base = {{1, 0, 0, 0, 0, 0, 1, 0}, {0.1, 0, 10, 1, 0, 0, 1, 0.1}, {1, 0.5, 1, 0.5, 0, 0, 1, -0.1}};
        for (int a = 0; a < mu + lam; a++) {
            Policy P; P.w.assign(NW, 0); P.reset = o_reset; P.curdef = o_curdef; P.Tw = o_Tw; P.hard = o_hard; P.stagger = o_stag; P.quota = o_quota; P.qfrom = g_qfrom; P.plan = o_plan; P.planf = o_planf; P.planN = g_planN; P.lpz = o_lpz;
            if (a < 3) { for (int q = 0; q < 8; q++) P.w[q] = base[a][q]; if (o_reset || o_pad || !o_lpz.empty()) P.w[8] = -1; if (!o_lpz.empty()) P.w[13] = 1.0; }   // 基础个体：从不提前大修
            else for (auto& x : P.w) x = N01(rng);
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
                    for (int q = 6; q < NW; q++) P.w[q] = P.w[q] + 0.3 * N01(rng) * (fabs(P.w[q]) + 0.2);
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
        printf("{\"inst\": \"%s\", \"mode\": \"%s\", \"H\": %lld, \"gp\": %d, \"reset\": %d, \"curdef\": %d, \"Tw\": %lld, \"hard\": %d, \"stagger\": %d, \"seed\": %llu, \"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, "
               "\"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, \"duties\": %lld, \"cyc_resets\": %lld, \"cyc_duties\": %lld, "
               "\"sec\": %.1f, \"hist\": %s, \"w\": %s, \"tree\": \"%s\"}\n",
               I.name.c_str(), mode.c_str(), H, gp ? 1 : 0, o_reset ? 1 : 0, o_curdef ? 1 : 0, o_Tw, o_hard ? 1 : 0, o_stag ? 1 : 0, seed, b.s.R, b.s.status.c_str(), vec_s(b.s.Nk).c_str(), vec_s(b.s.Mj).c_str(),
               v.gaps, v.cert, v.period, exact, v.resets, v.duties, v.cyc_resets, v.cyc_duties, sec, hs.c_str(), dvec_s(b.P.w).c_str(), gp ? G.str(b.P.tree, 0).c_str() : "");
        return 0;
    }
    return 1;
}
