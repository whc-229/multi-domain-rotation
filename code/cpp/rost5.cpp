// rost5（code/cpp/rost5.cpp，T3）：在 prior 的 rost4 上做三处扩展，其余逐字不变：
//   1) 每站 K_s ≥ 1：站 s 拆成 K_s 个岗位，每个岗位是一条首尾相接的班链（"每格恰好 K_s 台在岗"的等价形式）。K_s 全为 1 时与 rost4 完全相同。
//   2) rlp=<文件>：由松弛推导的派工规则 R_LP（见 RLP 结构体与 run() 里的 rlp 分支），可配 planf=<文件> 的 LP 提前大修计划；
//      文件里带 steady 标记时用稳态初值（各相位台数按权重分配、就绪时刻均匀错开）。
//   3) 新命令 rlp：从 LP 取整台数出发，沿用 rost4 的"补足 / 回删 + 证书"定台数循环，然后 Tverify 格复核。
//   用法补充：rost5 rlp <inst> <rlp文件> <H定台数> <Tcert> <Tverify> [planf=文件]
//            规则文件带 noscale 1 时（T3b 的 _u 变体）三种计划都不按台数缩放。
//            rost5 replay <inst> <排班文件> <Tmax>（T4：独立重放上界解的显式排班，检查每班合法、每格覆盖）
//   T6：定台数循环证书阶段的补台默认改为"试加每种资源各 1，取 Tcert 格内缺岗最少的"（同缺岗时取缺岗归因信号大的，再同取下标小的）；
//       当前台数的 Tcert 格缺岗连续 5 次补台没有刷新最低值时，判本起点失败（不再空转补到 60 次）。
//       旧做法（直接加归因信号最大的）用选项 certadd=why 保留；不带停止条件的试加用 certadd=trial60。rlp、es、size、rdd 命令都接受这两个选项。
//   T7：规则偏离解码 RDD（Policy::rdd）。窗口 [0, W) 内每一次交接（按处理顺序编号）是一个决策；
//       基因型给出若干 (决策序号, u)；该决策上把默认规则的排序里除第一名以外可行的机器人依次排成备选表
//       （允许时再加"本班后提前大修"一项），选第 ⌊u·|A|⌋ 项。全部槽关闭时与默认规则逐位相同。
//       新命令：rost5 size <inst> <H> <Tcert> <Tverify> <规则选项…>    给定规则（可带 rdd=文件）定台数 + 复核
//               rost5 rdd  <inst> <H> <Tcert> <Tverify> <gens> <mu> <lam> <seed> <rho> <qon> <规则选项…>   RDD 的 (μ+λ) 进化
//       规则选项：rlp=<文件> | w=<逗号分隔权重> [lpg=<文件>] [curdef] [planf=<文件>] [facloop]，以及 W=<格>、resetopt、certadd=why。
// 以下为 rost4 原注释。
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
        psite.clear(); for (int s = 0; s < ns; s++) for (int q = 0; q < S[s].K; q++) psite.push_back(s);   // 岗位 → 站
    }
    vector<int> psite;
};

// R_LP 规则数据（由 code/t3_make_rlp.py 生成，全部整数，保证与 Python 参考实现逐位一致）
struct RLP {
    bool on = false, steady = false; int H = 0;
    bool noscale = false;               // T3b（_u 变体）：值岗计划、相位计划、提前大修计划都不按 Nk/N_LP 缩放，直接用 LP 累计流
    vector<ll> NLPq;                    // LP 台数 × 1e6（四舍五入）
    vector<int> Q, CP;                  // 每站 lag 单位 1/Q_s；每型相位 lag 截断 ±CP_k（单位 1/L_k）
    vector<vector<int>> p, CQ;          // p[k][s]：t ≥ H 后每格计划增量（单位 1/Q_s）；CQ[k][s]：lag 截断
    vector<vector<vector<ll>>> F, FP;   // F[k][s][t]、FP[k][i][t]：t < H 的 LP 累计值岗流 × 1e6（D ≤ t）
    vector<vector<int>> wph;            // 稳态初值：每型每相位的权重（= 就绪窗口长度）
    vector<int> N0, M0;
};

// RDD 解码数据（T7）：dev[c] = 第 c 个窗口决策上的 u ∈ [0, 1)，−1 表示不偏离；maxpos = 最后一个偏离决策的序号（无偏离为 −1）
struct RDD { ll W = 0; int ndev = 0; vector<double> dev; int maxpos = -1; bool reset_opt = false; };

// ---------------- 策略：线性打分（6 维机器人 + 2 维工具）或 GP 树（机器人打分）+ 线性工具打分 ----------------
struct GPNode { int op; double c; int a, b; };   // op: 0..5 特征, 6 常数, 7 +, 8 -, 9 *, 10 保护除, 11 max, 12 min
struct Policy {
    vector<double> w;            // 至少 8 维；reset 时 13 维
    vector<GPNode> tree;         // 空则用线性
    bool reset = false, curdef = false, hard = false, stagger = false; ll Tw = 0;
    // T13：逐决策编码。pd_mode = 1：RK（原排班软件的随机键，可用机器人按编号）；2：D5（上一篇的 dec5，按出勤少、相位低、编号排序，带回溯修复）。
    // 窗口 [0, pdW) 内第 c 个决策取第 ⌊pdg[c]·n⌋ 个可行候选（n = 工具也够用的候选数；缺基因按 0）；pdf[c] ≥ 0 时改取第 pdf[c] 个（D5 修复的强制选择）。
    // 窗口之后按该解码器自己的排序取第一个（相当于基因取 0），不再修复。stop_dec / stop_wgap / wgap_from 只供修复的内部推演用。
    int pd_mode = 0; ll pdW = 0; const vector<double>* pdg = nullptr; const vector<int>* pdf = nullptr; bool repair = false;
    ll stop_dec = -1; bool stop_wgap = false; ll wgap_from = 0;
    ll plan_until = -1;          // T12：planW=，计划提前大修只在 t < plan_until 时执行（−1 不限；严格的"只在启动窗口内提前大修"）
    vector<vector<int>> quota; vector<vector<ll>> qfrom;
    vector<vector<vector<pair<ll,double>>>> planf; vector<double> planN;
    const RLP* rlp = nullptr;
    const RDD* rdd = nullptr;     // T7：规则偏离解码（nullptr 时与原规则完全相同）
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
             ll resets = 0, duties = 0, cyc_resets = -1, cyc_duties = -1; ll ndec_w = 0; ll last_reset = -1; int stop_ncand = -1; ll gap_dec = -1; };   // stop_ncand / gap_dec：T13 修复推演的停点信息   // last_reset：最后一次提前大修的时刻（−1 = 没有）   // ndec_w：窗口 [0, W) 内的决策数（RDD）
static inline pair<ull,ull> hsh(const vector<ll>& v) { ull a = 1469598103934665603ULL, b = 0x9E3779B97F4A7C15ULL;
  for (ll x : v) { ull y = (ull)x; a ^= y + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2); a *= 1099511628211ULL; b = (b ^ (y * 0xff51afd7ed558ccdULL)) * 0xc4ceb9fe1a85ec53ULL; b ^= b >> 33; }
  return {a, b}; }
struct PH { size_t operator()(const pair<ull,ull>& p) const { return p.first ^ (p.second * 31); } };

static inline int postu(const Inst& I, int k, int s, int j, double unin) {
    const Type& t = I.T[k];
    return (int)ceil((t.post_h + t.back_h + I.S[s].delta + unin + I.TT[j].maint) / 2.0 - 1e-9);
}

// T10 诊断：g_gaplog 非空时，每个缺岗写一行 "t 站 就绪候选数 归因（1 = 工具，0 = 机器人）"；g_sizelog 非空时，定台数循环每一步写一行。
// 两者只写文件，不改变任何决策。
static FILE* g_gaplog = nullptr; static FILE* g_sizelog = nullptr;
// gap_cap（T7 加速）：不判环的运行里，缺岗数一旦超过 gap_cap 就提前返回（此时只保证返回的 gaps > gap_cap）。
// 定台数循环里"试加取缺岗最少"只需要知道候选是否可能胜出，超过当前最好值的候选不可能被选中，所以选择结果与不设上限时逐位相同。
Res run_core(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll Tmax, bool certify, ll stop_at = -1, vector<ll>* snap = nullptr, ll gap_cap = LLONG_MAX) {
    int nk = I.nk, ns = I.ns, nt = I.nt;
    vector<int> rk, rph; vector<ll> rv;
    vector<ll> att;   // T13 D5：每台机器人的出勤（值岗）次数
    vector<int> grp; vector<char> firstc; vector<vector<int>> qleft = P.quota; auto pl = P.plan;
    if (!P.planf.empty()) {   // 按实际台数换算：累计流量 × Nk/N_LP 四舍五入
        pl.assign(nk, {});
        for (int k = 0; k < nk; k++) {
            pl[k].assign(I.T[k].L, {});
            double sc = (P.planN[k] > 1e-9) ? (double)Nk[k] / P.planN[k] : 0.0;
            if (P.rlp && P.rlp->on && P.rlp->noscale) sc = 1.0;   // _u：LP 累计流直接四舍五入成整数次
            for (int i = 0; i < I.T[k].L; i++) { double F = 0; int prev = 0;
                for (auto& e2 : P.planf[k][i]) { F += e2.second * sc; int c = (int)llround(F) - prev; if (c > 0) { pl[k][i].push_back({e2.first, c}); prev += c; } } }
        }
    }
    const RLP* RL = (P.rlp && P.rlp->on) ? P.rlp : nullptr;
    if (RL && RL->steady) {   // 稳态初值：相位 i 的台数 ∝ wph[k][i]（最大余数法），第 j 台就绪时刻 = j·wph/n_i
        for (int k = 0; k < nk; k++) {
            int L = I.T[k].L; ll W = 0; for (int i = 0; i < L; i++) W += RL->wph[k][i];
            vector<ll> n(L), rem(L); ll used = 0;
            for (int i = 0; i < L; i++) { n[i] = (ll)Nk[k] * RL->wph[k][i] / W; rem[i] = (ll)Nk[k] * RL->wph[k][i] % W; used += n[i]; }
            vector<int> ord(L); iota(ord.begin(), ord.end(), 0);
            stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return rem[a] > rem[b]; });
            for (ll e = 0; e < Nk[k] - used; e++) n[ord[e]]++;
            int c = 0;
            for (int i = 0; i < L; i++) for (ll j = 0; j < n[i]; j++, c++) { rk.push_back(k); rph.push_back(i); rv.push_back(j * RL->wph[k][i] / n[i]); grp.push_back(c % L); firstc.push_back(1); }
        }
    } else
    for (int k = 0; k < nk; k++) for (int c = 0; c < Nk[k]; c++) { rk.push_back(k); rph.push_back(0); rv.push_back(I.T[k].v0); grp.push_back(c % I.T[k].L); firstc.push_back(1); }
    int nr = rk.size();
    att.assign(nr, 0);
    bool use_tools = nt > 0;
    vector<vector<ll>> tf(nt);
    for (int j = 0; j < nt; j++) tf[j].assign(Mj[j], -1000000);
    int c0 = 0; for (int k = 0; k < nk; k++) c0 = max(c0, I.T[k].v0);
    const int np = I.psite.size();
    vector<ll> fr(np); for (int q = 0; q < np; q++) fr[q] = c0 + I.S[I.psite[q]].lead;
    // R_LP 的 lag 状态
    vector<vector<ll>> lag(nk, vector<ll>(ns, 0)), phl(nk);
    vector<vector<ll>> scaleQ(nk, vector<ll>(ns, 0));
    for (int k = 0; k < nk; k++) phl[k].assign(I.T[k].L, 0);
    auto planQ = [&](int k, int s, ll t) -> ll {   // 计划累计班数（单位 1/Q_s）
        if (t < 0) return 0;
        const RLP& R = *RL; ll Hh = R.H;
        auto tab = [&](ll tt) -> ll { if (R.F[k][s].empty()) return 0; if (R.noscale) return (ll)R.Q[s] * R.F[k][s][tt] / 1000000; if (R.NLPq[k] <= 0) return 0; return (ll)R.Q[s] * Nk[k] * R.F[k][s][tt] / R.NLPq[k]; };
        if (t < Hh) return tab(t);
        ll B = Hh > 0 ? tab(Hh - 1) : 0;
        return B + (ll)R.p[k][s] * (t - Hh + 1);
    };
    auto phplanQ = [&](int k, int i, ll t) -> ll {   // 相位计划（单位 1/L_k）；t ≥ H 后不随时间变（改由派工驱动）
        const RLP& R = *RL; ll Hh = R.H;
        if (t < 0 || Hh == 0) return 0;
        ll tt = min(t, Hh - 1);
        if (R.noscale) return (ll)I.T[k].L * R.FP[k][i][tt] / 1000000;
        if (R.NLPq[k] <= 0) return 0;
        return (ll)I.T[k].L * Nk[k] * R.FP[k][i][tt] / R.NLPq[k];
    };
    auto clip = [](ll x, ll c) { return x > c ? c : (x < -c ? -c : x); };
    ll tprev = -1;
    vector<vector<ll>> D(nk); vector<ll> tot(nk, 0);
    for (int k = 0; k < nk; k++) D[k].assign(I.T[k].L, 0);
    Res R; R.gaps = 0; R.cert = -1; R.period = -1; R.why_r.assign(nk, 0); R.why_t.assign(nt, 0);
    vector<vector<int>> cnt(nk); for (int k = 0; k < nk; k++) { cnt[k].assign(I.T[k].L, 0); cnt[k][0] = Nk[k]; }
    const bool use_tw = P.reset && P.Tw > 0;
    ll qmax = 0; for (auto& v : P.qfrom) for (ll x : v) qmax = max(qmax, x);
    for (auto& pk : P.plan) for (auto& pi : pk) for (auto& e2 : pi) qmax = max(qmax, e2.first);
    for (auto& pk : P.planf) for (auto& pi : pk) for (auto& e2 : pi) qmax = max(qmax, e2.first);   // 配额起始时刻之前策略依赖 t，不判环
    if (RL) qmax = max(qmax, (ll)RL->H);   // R_LP：t < H 时计划表依赖 t
    // 判环：状态第一次重复
    vector<ll> cur; unordered_map<pair<ull,ull>, ll, PH> seen; vector<ll> st_res, st_dut;
    vector<int> due; vector<pair<double,int>> cand; vector<tuple<double,ll,int,int>> opts;
    ll t = *min_element(fr.begin(), fr.end());
    const bool lin = P.tree.empty();
    const RDD* DV = P.rdd;
    // 机器人 r 此刻去 s 站是否有足够的空闲适配工具（与派工时的工具检查相同，但不改状态；RDD 构造备选表用）
    auto tool_ok = [&](int r, int s) -> bool {
        if (!use_tools) return true;
        int k = rk[r], e = I.eidx[k][s]; const Type& T = I.T[k]; ll need = t - T.pre[e]; int c = 0;
        for (int j : T.C) for (ll f : tf[j]) if (f <= need) c++;
        return c >= T.q;
    };
    // RDD：在按默认规则排好序的候选里取可行者 fe（第一名即默认选择），按 u 选备选项。
    // 返回 −1 = 不偏离；≥ 0 = 改派的机器人编号；−2 = 仍派默认机器人，但本班后提前大修
    auto pick_dev = [&](const vector<int>& order, int s, double u) -> int {
        vector<int> fe; for (int r : order) if (tool_ok(r, s)) fe.push_back(r);
        if (fe.empty()) return -1;
        int r0 = fe[0]; bool ro = DV->reset_opt && rph[r0] != I.T[rk[r0]].L - 1;
        int nA = (int)fe.size() - 1 + (ro ? 1 : 0);
        if (nA <= 0) return -1;
        int a = (int)floor(u * (double)nA); if (a >= nA) a = nA - 1;
        if (a < (int)fe.size() - 1) return fe[a + 1];
        return -2;
    };
    while (t < Tmax) {
        if (RL) {   // 计划增量计入 lag（截断）
            for (int k = 0; k < nk; k++) {
                for (int s : I.T[k].elig) lag[k][s] = clip(lag[k][s] + planQ(k, s, t) - planQ(k, s, tprev), RL->CQ[k][s]);
                for (int i = 0; i < I.T[k].L; i++) phl[k][i] = clip(phl[k][i] + phplanQ(k, i, t) - phplanQ(k, i, tprev), RL->CP[k]);
            }
            tprev = t;
        }
        if (certify) {
            cur.clear();
            for (int r = 0; r < nr; r++) { cur.push_back(rk[r]); cur.push_back(rph[r]); cur.push_back(max(rv[r] - t, (ll)-CL)); if (P.stagger) cur.push_back(firstc[r]); }
            for (auto& qk : qleft) for (int q : qk) cur.push_back(q);
            for (auto& pk : pl) for (auto& pi : pk) for (auto& e2 : pi) cur.push_back(e2.second);
            for (int q = 0; q < np; q++) cur.push_back(fr[q] - t);
            if (RL) { for (int k = 0; k < nk; k++) { for (int s : I.T[k].elig) cur.push_back(lag[k][s]); for (ll x : phl[k]) cur.push_back(x); } }
            if (!P.curdef) for (int k = 0; k < nk; k++) for (int i = 0; i < I.T[k].L; i++) cur.push_back(tot[k] - (ll)I.T[k].L * D[k][i]);
            if (P.pd_mode == 2 && nr > 0) { ll mn = *min_element(att.begin(), att.end()); for (int r = 0; r < nr; r++) cur.push_back(min(att[r] - mn, (ll)CL)); }
            for (int j = 0; j < nt; j++) {
                size_t st0 = cur.size();
                for (ll f : tf[j]) cur.push_back(max(f - t, (ll)-CL));
                sort(cur.begin() + st0, cur.end());
            }
            if (stop_at >= 0) { if (t >= stop_at) { *snap = cur; R.t_end = t; return R; } }
            else if ((!use_tw || t >= P.Tw) && t >= qmax && (!DV || t >= DV->W || R.ndec_w > DV->maxpos) && (!P.pd_mode || t >= P.pdW)) {   // RDD：最后一个偏离决策之后才判环
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
        for (int q = 0; q < np; q++) if (fr[q] <= t) due.push_back(q);
        stable_sort(due.begin(), due.end(), [&](int a, int b) { return I.nelig[I.psite[a]] < I.nelig[I.psite[b]]; });
        for (int pq : due) {
            const int s = I.psite[pq];
            double dv = -1; bool force_reset = false;      // RDD：本决策的 u（−1 不偏离）、是否强制本班后提前大修
            if (DV && t < DV->W) { ll c = R.ndec_w++; if (c < DV->ndev) dv = DV->dev[c]; }
            ll pdc = -1; if (P.pd_mode && t < P.pdW) pdc = R.ndec_w++;   // T13：逐决策编码的窗口决策序号
            if (RL) {   // ---------- R_LP 派工 ----------
                vector<tuple<ll,int,ll,ll,int>> rc;   // (−lag_{k,s}, k, −相位 lag, 截断就绪, 编号)
                for (int r = 0; r < nr; r++) {
                    int k = rk[r], e = I.eidx[k][s];
                    if (e < 0 || rv[r] + I.S[s].lead > t) continue;
                    rc.emplace_back(-lag[k][s], k, -phl[k][rph[r]], max(rv[r] - t, (ll)-CL), r);
                }
                sort(rc.begin(), rc.end());
                if (dv >= 0) {
                    vector<int> order; for (auto& cr : rc) order.push_back(get<4>(cr));
                    int pk = pick_dev(order, s, dv);
                    if (pk == -2) force_reset = true;
                    else if (pk >= 0) { auto it = find_if(rc.begin(), rc.end(), [&](const tuple<ll,int,ll,ll,int>& x) { return get<4>(x) == pk; }); rotate(rc.begin(), it, it + 1); }
                }
                bool done = false;
                for (auto& cr : rc) {
                    int r = get<4>(cr), k = rk[r], i = rph[r], e = I.eidx[k][s]; const Type& T = I.T[k];
                    if (use_tools) {
                        ll need = t - T.pre[e];
                        vector<tuple<double,ll,int,int>> op;
                        for (int j : T.C) for (int idx = 0; idx < (int)tf[j].size(); idx++) {
                            ll f = tf[j][idx];
                            if (f <= need) op.emplace_back(I.TT[j].maint, max(f - t, (ll)-CL), j, idx);
                        }
                        if ((int)op.size() < T.q) continue;
                        partial_sort(op.begin(), op.begin() + T.q, op.end());
                        double unin = 0; for (int c = 0; c < T.q; c++) unin += I.TT[get<2>(op[c])].unin;
                        for (int c = 0; c < T.q; c++) { int j = get<2>(op[c]), idx = get<3>(op[c]); tf[j][idx] = t + T.d + postu(I, k, s, j, unin); }
                    }
                    int ld = T.step[i][e].first, nx = T.step[i][e].second, nph = (i + 1) % T.L;
                    bool planned = false;
                    if (!pl.empty() && i != T.L - 1 && (P.plan_until < 0 || t < P.plan_until)) {
                        for (auto& e2 : pl[k][i]) { if (e2.first > t) break; if (e2.second > 0) { e2.second--; planned = true; break; } }
                    }
                    if (planned) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; }
                    if (force_reset && nph != 0) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; }
                    if (nph == 0 && i != T.L - 1) R.last_reset = t;
                    rv[r] = t - ld + nx; rph[r] = nph;
                    cnt[k][i]--; cnt[k][nph]++; R.duties++;
                    D[k][i]++; tot[k]++;
                    lag[k][s] = clip(lag[k][s] - RL->Q[s], RL->CQ[k][s]);
                    phl[k][i] -= T.L;
                    if (t >= RL->H) for (int i2 = 0; i2 < T.L; i2++) phl[k][i2] += 1;
                    for (int i2 = 0; i2 < T.L; i2++) phl[k][i2] = clip(phl[k][i2], RL->CP[k]);
                    fr[pq] = t + T.d; done = true; break;
                }
                if (!done) {
                    R.gaps++;
                    if (g_gaplog) fprintf(g_gaplog, "%lld %d %d %d\n", t, s, (int)rc.size(), rc.empty() ? 0 : 1);
                    if (!rc.empty()) { for (int j : I.T[rk[get<4>(rc[0])]].C) R.why_t[j]++; }
                    else { for (int k = 0; k < nk; k++) if (I.eidx[k][s] >= 0) R.why_r[k]++; }
                    fr[pq] = t + 1;
                    if (!certify && R.gaps > gap_cap) { R.t_end = t; return R; }
                }
                continue;
            }
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
            if (dv >= 0) {
                vector<int> order; for (auto& cr : cand) order.push_back(cr.second);
                int pk = pick_dev(order, s, dv);
                if (pk == -2) force_reset = true;
                else if (pk >= 0) { auto it = find_if(cand.begin(), cand.end(), [&](const pair<double,int>& x) { return x.second == pk; }); rotate(cand.begin(), it, it + 1); }
            }
            if (P.pd_mode) {   // T13：按解码器自己的顺序重排候选，窗口内按基因（或修复的强制选择）取一个可行候选
                if (P.pd_mode == 1) sort(cand.begin(), cand.end(), [](const pair<double,int>& a, const pair<double,int>& b) { return a.second < b.second; });
                else {
                    ll mn = *min_element(att.begin(), att.end());
                    sort(cand.begin(), cand.end(), [&](const pair<double,int>& a, const pair<double,int>& b) {
                        ll xa = min(att[a.second] - mn, (ll)CL), xb = min(att[b.second] - mn, (ll)CL);
                        if (xa != xb) return xa < xb; if (rph[a.second] != rph[b.second]) return rph[a.second] < rph[b.second]; return a.second < b.second; });
                }
                if (pdc >= 0) {
                    vector<int> fe; for (auto& cr : cand) if (tool_ok(cr.second, s)) fe.push_back(cr.second);
                    if (pdc == P.stop_dec) { R.stop_ncand = (int)fe.size(); R.t_end = t; return R; }
                    if (fe.empty()) { if (P.stop_wgap && pdc >= P.wgap_from) { R.gap_dec = pdc; R.t_end = t; return R; } }
                    else {
                        int n = (int)fe.size(), idx;
                        int fc = (P.pdf && pdc < (ll)P.pdf->size()) ? (*P.pdf)[pdc] : -1;
                        if (fc >= 0) idx = min(fc, n - 1);
                        else { double x = (P.pdg && pdc < (ll)P.pdg->size()) ? (*P.pdg)[pdc] : 0.0; idx = (int)floor(x * (double)n); if (idx >= n) idx = n - 1; if (idx < 0) idx = 0; }
                        int pk = fe[idx]; auto it = find_if(cand.begin(), cand.end(), [&](const pair<double,int>& x) { return x.second == pk; }); rotate(cand.begin(), it, it + 1);
                    }
                }
            }
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
                if (!pl.empty() && i != T.L - 1 && (P.plan_until < 0 || t < P.plan_until)) {
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
                if (force_reset && nph != 0) { nx = T.step[T.L - 1][e].second; nph = 0; R.resets++; }
                if (nph == 0 && i != T.L - 1) R.last_reset = t;
                rv[r] = t - ld + nx; rph[r] = nph; att[r]++;
                cnt[k][i]--; cnt[k][nph]++; R.duties++;
                D[k][i]++; tot[k]++;
                fr[pq] = t + T.d; done = true; break;
            }
            if (!done) {
                R.gaps++;
                if (g_gaplog) fprintf(g_gaplog, "%lld %d %d %d\n", t, s, (int)cand.size(), cand.empty() ? 0 : 1);
                if (!cand.empty()) { for (int j : I.T[rk[cand[0].second]].C) R.why_t[j]++; }
                else { for (int k = 0; k < nk; k++) if (I.eidx[k][s] >= 0) R.why_r[k]++; }
                fr[pq] = t + 1;
                if (!certify && R.gaps > gap_cap) { R.t_end = t; return R; }
            }
        }
        t = *min_element(fr.begin(), fr.end());
    }
    R.t_end = t; return R;
}

// T13 D5 的回溯修复（移植 dec5 的 robot_arrangement_bt，BT_WINDOW = 12 个决策）：只在窗口 [0, pdW) 内。
// 按顺序找窗口内的缺岗决策 c（没有可行候选）；从 c−1 往回到 c−12，对每个决策 tt 依次试别的候选 alt（跳过 tt 上已有的强制选择），
// tt 之后的强制选择清空、其余按基因重新推进到 c：若 [tt, c) 内没有缺岗且到 c 时有可行候选，就采用（强制 tt 取 alt）。都不行时 c 照常缺岗。
// 返回强制选择表（下标 = 窗口决策序号，−1 = 不强制）。结果只取决于台数与基因，所以解码仍是确定性的。
vector<int> d5_repair(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P) {
    const int BT = 12; vector<int> f; Policy Q = P; Q.repair = false; ll from = 0;
    while (true) {
        Q.pdf = &f; Q.stop_dec = -1; Q.stop_wgap = true; Q.wgap_from = from;
        Res r = run_core(I, Nk, Mj, Q, P.pdW, false);
        if (r.gap_dec < 0) break;
        ll c = r.gap_dec; bool fixed = false;
        for (ll tt = c - 1; tt >= max((ll)0, c - BT) && !fixed; tt--) {
            vector<int> f2 = f; if ((ll)f2.size() < c + 1) f2.resize(c + 1, -1);
            for (ll j = tt; j < (ll)f2.size(); j++) f2[j] = -1;
            Q.pdf = &f2; Q.stop_wgap = false; Q.stop_dec = tt;
            Res r2 = run_core(I, Nk, Mj, Q, P.pdW, false);
            if (r2.stop_ncand <= 0) continue;
            int n = r2.stop_ncand, cur = (tt < (ll)f.size()) ? f[tt] : -1;
            for (int alt = 0; alt < n && !fixed; alt++) {
                if (alt == cur) continue;
                f2[tt] = alt; Q.stop_dec = c;
                Res r3 = run_core(I, Nk, Mj, Q, P.pdW, false);
                if (r3.stop_ncand > 0 && r3.gaps == r2.gaps) { f = f2; fixed = true; }
            }
        }
        from = c + 1;
    }
    return f;
}
// 仿真入口：D5 带修复时先算强制选择表，再仿真；其余情况就是 run_core。
Res run(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll Tmax, bool certify, ll stop_at = -1, vector<ll>* snap = nullptr, ll gap_cap = LLONG_MAX) {
    if (P.pd_mode == 2 && P.repair) { vector<int> f = d5_repair(I, Nk, Mj, P); Policy Q = P; Q.repair = false; Q.pdf = &f; return run_core(I, Nk, Mj, Q, Tmax, certify, stop_at, snap, gap_cap); }
    return run_core(I, Nk, Mj, P, Tmax, certify, stop_at, snap, gap_cap);
}

bool exact_check(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll enter, ll period) {
    vector<ll> a, b;
    Res r1 = run(I, Nk, Mj, P, LLONG_MAX / 4, true, enter, &a);
    Res r2 = run(I, Nk, Mj, P, LLONG_MAX / 4, true, enter + period, &b);
    return r1.t_end == enter && r2.t_end == enter + period && a == b;
}

struct SizeRes { int R; vector<int> Nk, Mj; string status; ll period; };
string vec_s(const vector<int>& v);

vector<int> g_N0, g_M0;   // LP 引导时的定台数起点
// 证书阶段补台：1 = 试加每种资源取 Tcert 格内缺岗最少；当前台数的 Tcert 格缺岗连续 CERT_PATIENCE 次补台都没有刷新最低值时判本起点失败（T6 起的默认）；
//               2 = 试加但不设停止条件，最多补 max_add 次（certadd=trial60，T6 开发中间版本，只为复现用）；0 = 直接加归因信号最大的（旧做法，certadd=why）
int g_certadd = 1; const int CERT_PATIENCE = 5;
SizeRes size_res(const Inst& I, const Policy& P, ll H, const string& mode, ll Tcert, int max_add = 60) {
    vector<int> Nk = g_N0.empty() ? I.N0 : g_N0, Mj = g_M0.empty() ? I.M0 : g_M0; int nk = I.nk, nt = I.nt, nres = nk + nt;
    auto bump = [&](int idx, int d) { if (idx < nk) Nk[idx] += d; else Mj[idx - nk] += d; };
    auto tot = [&]() { int s = 0; for (int x : Nk) s += x; for (int x : Mj) s += x; return s; };
    Res r = run(I, Nk, Mj, P, H, false); int added = 0;
    while (r.gaps > 0 && added < max_add) {
        vector<ll> sig = r.why_r; sig.insert(sig.end(), r.why_t.begin(), r.why_t.end());
        ll bg = LLONG_MAX, bsig = 0; int bi = -1;
        for (int idx = 0; idx < nres; idx++) {
            bump(idx, 1); ll g = run(I, Nk, Mj, P, H, false, -1, nullptr, bg).gaps; bump(idx, -1);   // 超过当前最好值即停（不影响选择）
            if (bi < 0 || g < bg || (g == bg && sig[idx] > bsig)) { bg = g; bi = idx; bsig = sig[idx]; }
        }
        bump(bi, 1); added++; r = run(I, Nk, Mj, P, H, false);
        if (g_sizelog) fprintf(g_sizelog, "fill add=%d gaps=%lld Nk=%s Mj=%s\n", bi, r.gaps, vec_s(Nk).c_str(), vec_s(Mj).c_str());
    }
    if (r.gaps > 0) { if (g_sizelog) fprintf(g_sizelog, "fail_fill\n"); return {-1, Nk, Mj, "fail", -1}; }
    if (g_sizelog) fprintf(g_sizelog, "filled Nk=%s Mj=%s\n", vec_s(Nk).c_str(), vec_s(Mj).c_str());
    bool improved = true;
    while (improved) {
        improved = false;
        for (int idx = 0; idx < nres; idx++) {
            int& v = idx < nk ? Nk[idx] : Mj[idx - nk];
            if (v == 0) continue;
            v--;
            if (run(I, Nk, Mj, P, H, false, -1, nullptr, 0).gaps == 0) { improved = true; break; }   // 出现第 1 个缺岗即停
            v++;
        }
    }
    if (g_sizelog) fprintf(g_sizelog, "pruned Nk=%s Mj=%s\n", vec_s(Nk).c_str(), vec_s(Mj).c_str());
    if (mode == "finite") return {tot(), Nk, Mj, "finite", -1};
    bool strict = (mode == "certs");
    ll best_g0 = LLONG_MAX; int stall = 0;      // g_certadd = 1 的停止条件：缺岗最低值与连续未刷新的次数
    for (int a = 0; a < max_add; a++) {
        Res c = run(I, Nk, Mj, P, Tcert, true);
        if (strict && c.cert == -1 && c.gaps == 0) c = run(I, Nk, Mj, P, 4 * Tcert, true);   // 判不出循环：加长到 4 倍再判
        if (g_sizelog) fprintf(g_sizelog, "cert a=%d gaps=%lld cert=%d period=%lld t_end=%lld Nk=%s Mj=%s\n", a, c.gaps, c.cert, c.period, c.t_end, vec_s(Nk).c_str(), vec_s(Mj).c_str());
        if (c.cert == 1) return {tot(), Nk, Mj, "cert", c.period};
        if (c.cert == -1 && c.gaps == 0) return {tot(), Nk, Mj, "uncert", -1};
        vector<ll> sig = c.why_r; sig.insert(sig.end(), c.why_t.begin(), c.why_t.end());
        if (g_certadd == 0) { int bi = max_element(sig.begin(), sig.end()) - sig.begin(); bump(bi, 1); continue; }
        // 试加：每种资源各加 1，跑满 Tcert 格（不判环）数缺岗，取最少的；同缺岗取归因信号大的，再同取下标小的（与补足阶段相同）。
        // 默认（g_certadd = 1）：记当前台数的 Tcert 格缺岗 g0；连续 CERT_PATIENCE 次补台 g0 都没有低于此前最低值，判本起点失败（继续补多半是空转）。
        if (g_certadd == 1) {
            ll g0 = run(I, Nk, Mj, P, Tcert, false).gaps;
            if (g_sizelog) fprintf(g_sizelog, "cert_g0 a=%d g0=%lld\n", a, g0);
            if (g0 < best_g0) { best_g0 = g0; stall = 0; } else if (++stall >= CERT_PATIENCE) { if (g_sizelog) fprintf(g_sizelog, "fail_cert_stall\n"); return {-1, Nk, Mj, "fail", -1}; }
        }
        ll bg = LLONG_MAX, bsig = 0; int bi = -1;
        for (int idx = 0; idx < nres; idx++) {
            bump(idx, 1); ll g = run(I, Nk, Mj, P, Tcert, false, -1, nullptr, bg).gaps; bump(idx, -1);   // 超过当前最好值即停（不影响选择）
            if (bi < 0 || g < bg || (g == bg && sig[idx] > bsig)) { bg = g; bi = idx; bsig = sig[idx]; }
        }
        bump(bi, 1);
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
// R_LP 文件：见 code/t3_make_rlp.py 的 write_rlp()
RLP parse_rlp(const Inst& I, const string& fn) {
    RLP R; R.on = true; ifstream f(fn); if (!f) { fprintf(stderr, "cannot open %s\n", fn.c_str()); exit(1); }
    string tag; int nk = I.nk, ns = I.ns;
    R.p.assign(nk, vector<int>(ns, 0)); R.CQ.assign(nk, vector<int>(ns, 0)); R.F.assign(nk, vector<vector<ll>>(ns));
    R.FP.assign(nk, {}); for (int k = 0; k < nk; k++) R.FP[k].assign(I.T[k].L, {});
    R.wph.assign(nk, {});
    while (f >> tag) {
        if (tag == "H") f >> R.H;
        else if (tag == "steady") { int x; f >> x; R.steady = x; }
        else if (tag == "noscale") { int x; f >> x; R.noscale = x; }
        else if (tag == "NLPq") { R.NLPq.resize(nk); for (auto& x : R.NLPq) f >> x; }
        else if (tag == "N0") { R.N0.resize(nk); for (auto& x : R.N0) f >> x; }
        else if (tag == "M0") { R.M0.resize(I.nt); for (auto& x : R.M0) f >> x; }
        else if (tag == "Q") { R.Q.resize(ns); for (auto& x : R.Q) f >> x; }
        else if (tag == "CP") { R.CP.resize(nk); for (auto& x : R.CP) f >> x; }
        else if (tag == "p") { int k, s; f >> k >> s; f >> R.p[k][s] >> R.CQ[k][s]; }
        else if (tag == "W") { int k; f >> k; R.wph[k].resize(I.T[k].L); for (auto& x : R.wph[k]) f >> x; }
        else if (tag == "F") { int k, s; f >> k >> s; R.F[k][s].resize(R.H); for (auto& x : R.F[k][s]) f >> x; }
        else if (tag == "FP") { int k, i; f >> k >> i; R.FP[k][i].resize(R.H); for (auto& x : R.FP[k][i]) f >> x; }
        else { fprintf(stderr, "bad rlp tag %s\n", tag.c_str()); exit(1); }
    }
    for (int k = 0; k < nk; k++) for (int i = 0; i < I.T[k].L; i++) if ((int)R.FP[k][i].size() != R.H) R.FP[k][i].assign(R.H, 0);
    return R;
}


// ---------------- T7：RDD 基因型文件、规则选项、定台数 ----------------
// 基因型文件：行 "W 格"、"ndev n"、"qon x"、"resetopt 0/1"，其余每行 "slot p a"（按槽的顺序）
struct Geno { ll W = 0; int ndev = 0; double qon = 0.4; bool reset_opt = false; vector<double> g; };   // g = p0 a0 p1 a1 …
Geno read_geno(const string& fn) {
    Geno G; ifstream f(fn); if (!f) { fprintf(stderr, "cannot open %s\n", fn.c_str()); exit(1); }
    string tag; while (f >> tag) {
        if (tag == "W") f >> G.W; else if (tag == "ndev") f >> G.ndev; else if (tag == "qon") f >> G.qon;
        else if (tag == "resetopt") { int x; f >> x; G.reset_opt = x; }
        else if (tag == "slot") { double p, a; f >> p >> a; G.g.push_back(p); G.g.push_back(a); }
        else { fprintf(stderr, "bad geno tag %s\n", tag.c_str()); exit(1); }
    }
    return G;
}
// 解码：槽 (p, a) 在 a ≥ 1 − q_on 时激活，作用于第 ⌊p·N_dev⌋ 个窗口决策，u = (a − 1 + q_on)/q_on；同一决策上后一个槽生效
RDD decode_geno(const vector<double>& g, ll W, int ndev, double qon, bool reset_opt) {
    RDD D; D.W = W; D.ndev = ndev; D.reset_opt = reset_opt; D.dev.assign(max(ndev, 0), -1.0); D.maxpos = -1;
    const double thr = 1.0 - qon;
    for (size_t a = 0; a + 1 < g.size(); a += 2) {
        double p = g[a], x = g[a + 1];
        if (ndev <= 0 || x < thr) continue;
        int pos = (int)floor(p * (double)ndev); if (pos >= ndev) pos = ndev - 1;
        D.dev[pos] = (x - 1.0 + qon) / qon; D.maxpos = max(D.maxpos, pos);
    }
    return D;
}
// 规则选项：rlp=文件（R_LP）| w=逗号分隔权重 [lpg=文件] [curdef] [planf=文件]；另有 W=格、resetopt、certadd=why、rdd=基因型文件
struct Rule { Policy P; bool is_rlp = false, facloop = false; ll W = 0; bool reset_opt = false; string geno; };
extern ll g_hsfb;
// T11 定台数下降（rddfd）的选项：fdsel = 选减哪个资源时用的基因型（best = 最近一个判定可行的个体，base = 默认规则）；
// stall = 当前目标台数下连续多少代没有判定可行的解就换一个资源去减；fddump = 抽样记录（台数, 基因型, T_eval 缺岗）的文件，供核对
string g_fdsel = "best", g_fddump; int g_stall = 30; double g_swap = 0.5;   // swap：子代以此概率做一次台数互换（某资源 −1、另一资源 +1，总数不变）
static RLP g_RR;
Rule parse_rule(const Inst& I, int argc, char** argv, int a0) {
    Rule U; U.P.w.assign(13, 0);
    for (int a = a0; a < argc; a++) {
        string s = argv[a];
        if (s.rfind("rlp=", 0) == 0) { g_RR = parse_rlp(I, s.substr(4)); U.P.rlp = &g_RR; U.P.curdef = true; U.is_rlp = true; }
        else if (s.rfind("w=", 0) == 0) { U.P.w.clear(); stringstream ss(s.substr(2)); string x; while (getline(ss, x, ',')) U.P.w.push_back(atof(x.c_str())); }
        else if (s.rfind("lpg=", 0) == 0) U.P.lpz = parse_lpz(I, s.substr(4), g_N0, g_M0);
        else if (s == "curdef") U.P.curdef = true;
        else if (s.rfind("planf=", 0) == 0) { U.P.planf = parse_planf(I, s.substr(6)); U.P.planN = g_planN; U.P.curdef = true; }
        else if (s.rfind("W=", 0) == 0) U.W = atoll(s.c_str() + 2);
        else if (s == "resetopt") U.reset_opt = true;
        else if (s == "facloop") U.facloop = true;
        else if (s == "certadd=why") g_certadd = 0;
        else if (s == "certadd=trial60") g_certadd = 2;
        else if (s.rfind("rdd=", 0) == 0) U.geno = s.substr(4);
        else if (s.rfind("fdsel=", 0) == 0) g_fdsel = s.substr(6);
        else if (s.rfind("stall=", 0) == 0) g_stall = atoi(s.c_str() + 6);
        else if (s.rfind("swap=", 0) == 0) g_swap = atof(s.c_str() + 5);
        else if (s.rfind("fddump=", 0) == 0) g_fddump = s.substr(7);
        else if (s.rfind("sizelog=", 0) == 0) g_sizelog = fopen(s.substr(8).c_str(), "w");
        else if (s.rfind("gaplog=", 0) == 0) g_gaplog = fopen(s.substr(7).c_str(), "w");
        else if (s.rfind("hsfb=", 0) == 0) g_hsfb = atoll(s.c_str() + 5);
        else if (s.rfind("planW=", 0) == 0) U.P.plan_until = atoll(s.c_str() + 6);
        else if (s == "pd=rk") { U.P.pd_mode = 1; U.P.curdef = true; }
        else if (s == "pd=d5") { U.P.pd_mode = 2; U.P.repair = true; U.P.curdef = true; }
        else if (s == "pd=d5norepair") { U.P.pd_mode = 2; U.P.curdef = true; }
        else if (s.rfind("pdg=", 0) == 0) { static vector<double> G; G.clear(); stringstream ss(s.substr(4)); string x; while (getline(ss, x, ',')) G.push_back(atof(x.c_str())); U.P.pdg = &G; }
        else { fprintf(stderr, "bad rule option %s\n", s.c_str()); exit(1); }
    }
    while (U.P.w.size() < (U.P.lpz.empty() ? 13u : 14u)) U.P.w.push_back(0);
    if (U.P.pd_mode) U.P.pdW = U.W;
    return U;
}
// 按规则定台数：R_LP 与 rlp 命令相同（起点 LP 取整 × 1.00 / 1.15 / 1.30 / 1.50，取第一个成功的）；线性规则与 es 命令相同（起点 g_N0 或实例默认）。
// facloop：线性规则也用放大起点循环（起点基数为 g_N0 / g_M0，未设时用实例默认 N0 / M0）。调用后 g_N0 / g_M0 可能被改写，由调用方恢复。
// T10：hsfb = 补足阶段时域的后备值（0 = 不用）。R_LP 的 4 个起点在 Hs 格下都失败时，用 Hs = g_hsfb 把 4 个起点再跑一遍
// （m = 20 上补足时域 365 格短于规则依赖时刻的部分 H = 4·S_max，计划切换处的缺岗补足阶段看不到）。原来能定出台数的实例结果不变。
ll g_hsfb = 0; int g_hsfb_used = 0;
SizeRes size_rule_hs(const Inst& I, const Policy& P, bool is_rlp, ll Hs, ll Tcert, int* fac_used, bool facloop);
SizeRes size_rule(const Inst& I, const Policy& P, bool is_rlp, ll Hs, ll Tcert, int* fac_used = nullptr, bool facloop = false) {
    g_hsfb_used = 0;
    vector<int> n0 = g_N0, m0 = g_M0;
    SizeRes sr = size_rule_hs(I, P, is_rlp, Hs, Tcert, fac_used, facloop);
    if (sr.R < 0 && is_rlp && g_hsfb > 0 && g_hsfb != Hs) { g_N0 = n0; g_M0 = m0; g_hsfb_used = 1; sr = size_rule_hs(I, P, is_rlp, g_hsfb, Tcert, fac_used, facloop); }
    return sr;
}
SizeRes size_rule_hs(const Inst& I, const Policy& P, bool is_rlp, ll Hs, ll Tcert, int* fac_used, bool facloop) {
    if (!is_rlp && !facloop) { if (fac_used) *fac_used = 100; return size_res(I, P, Hs, "certs", Tcert); }
    vector<int> bN = is_rlp ? P.rlp->N0 : (g_N0.empty() ? I.N0 : g_N0), bM = is_rlp ? P.rlp->M0 : (g_M0.empty() ? I.M0 : g_M0);
    const int facs[4] = {100, 115, 130, 150}; SizeRes sr;
    for (int fi = 0; fi < 4; fi++) {
        g_N0.assign(I.nk, 0); g_M0.assign(I.nt, 0);
        for (int k = 0; k < I.nk; k++) g_N0[k] = (bN[k] * facs[fi] + 99) / 100;
        for (int j = 0; j < I.nt; j++) g_M0[j] = (bM[j] * facs[fi] + 99) / 100;
        sr = size_res(I, P, Hs, "certs", Tcert);
        if (sr.R >= 0) { if (fac_used) *fac_used = facs[fi]; return sr; }
    }
    if (fac_used) *fac_used = -1;
    return sr;
}
// 同资源数时的次级比较：从最终台数里任减 1 个资源后，Hs 格内缺岗数的最小值（越小越接近再省一个资源）
ll near_miss(const Inst& I, const Policy& P, const SizeRes& sr, ll Hs) {
    vector<int> Nk = sr.Nk, Mj = sr.Mj; ll best = LLONG_MAX;
    for (int idx = 0; idx < I.nk + I.nt; idx++) {
        int& v = idx < I.nk ? Nk[idx] : Mj[idx - I.nk];
        if (v == 0) continue;
        v--; best = min(best, run(I, Nk, Mj, P, Hs, false, -1, nullptr, best).gaps); v++;
    }
    return best == LLONG_MAX ? 0 : best;
}

// ---------------- T11：固定台数的完整判定（与 size 的判定部分相同） ----------------
// 6 万格判环（Tcert）→ 判不出且无缺岗时 24 万格（4·Tcert）→ Tverify 格复核（判环，判出时做 exact_check）。
// 与 size：size_res 证书阶段第一步（不补台）的判环，加 size 命令的复核，逐步相同。
// 可行 = 复核严格证明（verify_cert = 1 且 exact = 1）或复核 0 缺岗；不可行时 gaps 取判定中的缺岗数（判环阶段的，或复核的）。
struct Judge { bool feas = false; string status; ll gaps = 0; int cert = -1; ll vgaps = -1; int vcert = -2; ll vperiod = -1; int exact = -1; ll resets = 0, last_reset = -1; };
Judge judge_fleet(const Inst& I, const vector<int>& Nk, const vector<int>& Mj, const Policy& P, ll Tcert, ll Tver) {
    Judge J; Res c = run(I, Nk, Mj, P, Tcert, true);
    if (c.cert == -1 && c.gaps == 0) c = run(I, Nk, Mj, P, 4 * Tcert, true);
    J.cert = c.cert; J.gaps = c.gaps;
    if (c.cert == 1) J.status = "cert"; else if (c.cert == -1 && c.gaps == 0) J.status = "uncert";
    else { J.status = "fail"; return J; }
    Res v = run(I, Nk, Mj, P, Tver, true);
    J.vgaps = v.gaps; J.vcert = v.cert; J.vperiod = v.period; J.resets = v.resets; J.last_reset = v.last_reset;
    J.exact = (v.cert == 1) ? (exact_check(I, Nk, Mj, P, v.enter, v.period) ? 1 : 0) : -1;
    J.feas = (v.cert == 1 && J.exact == 1) || v.gaps == 0;
    if (!J.feas) J.gaps = max(v.gaps, (ll)1);
    return J;
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
            else { if (s.rfind("rlp=", 0) == 0) { static RLP RR; RR = parse_rlp(I, s.substr(4)); P.rlp = &RR; P.curdef = true; continue; }
                   if (s == "reset") P.reset = true; else if (s == "curdef") P.curdef = true; else if (s == "hard") P.hard = true; else if (s == "stagger") { P.stagger = true; P.curdef = true; } else if (s.rfind("quota=", 0) == 0) { P.curdef = true; P.quota = parse_quota(I, s.substr(6)); P.qfrom = g_qfrom; } else if (s.rfind("plan=", 0) == 0) { P.curdef = true; P.plan = parse_plan(I, s.substr(5)); } else if (s.rfind("planf=", 0) == 0) { P.curdef = true; P.planf = parse_planf(I, s.substr(6)); P.planN = g_planN; } else if (s.rfind("Tw=", 0) == 0) P.Tw = atoll(s.c_str() + 3);
                   else if (s.rfind("lpg=", 0) == 0) { static vector<int> n0, m0; P.lpz = parse_lpz(I, s.substr(4), n0, m0); }
                   else if (s.rfind("gaplog=", 0) == 0) g_gaplog = fopen(s.substr(7).c_str(), "w");
                   else if (s.rfind("planW=", 0) == 0) P.plan_until = atoll(s.c_str() + 6);
                   else if (s == "pd=rk") { P.pd_mode = 1; P.curdef = true; } else if (s == "pd=d5") { P.pd_mode = 2; P.repair = true; P.curdef = true; }
                   else if (s.rfind("pdW=", 0) == 0) P.pdW = atoll(s.c_str() + 4);
                   else if (s.rfind("pdg=", 0) == 0) { static vector<double> G; G.clear(); stringstream ss(s.substr(4)); string x; while (getline(ss, x, ',')) G.push_back(atof(x.c_str())); P.pdg = &G; }
                   else if (s.rfind("rdd=", 0) == 0) { static RDD DD; Geno G = read_geno(s.substr(4)); DD = decode_geno(G.g, G.W, G.ndev, G.qon, G.reset_opt); P.rdd = &DD; } }
        }
        if (P.reset) P.curdef = true;
        while (w.size() < (P.lpz.empty() ? 13u : 14u)) w.push_back(0);
        P.w = w;
        Res r = run(I, Nk, Mj, P, Tmax, cert);
        int exact = (r.cert == 1) ? (exact_check(I, Nk, Mj, P, r.enter, r.period) ? 1 : 0) : -1;
        printf("{\"gaps\": %lld, \"cert\": %d, \"period\": %lld, \"t_end\": %lld, \"exact\": %d, \"resets\": %lld, \"duties\": %lld, \"cyc_resets\": %lld, \"cyc_duties\": %lld, \"ndec_w\": %lld, \"last_reset\": %lld}\n",
               r.gaps, r.cert, r.period, r.t_end, exact, r.resets, r.duties, r.cyc_resets, r.cyc_duties, r.ndec_w, r.last_reset);
        return 0;
    }
    if (cmd == "replay") {   // T4：独立重放显式排班。文件：首行 "N …"、次行 "M …"，其余每行 "t s 机器人 提前大修 q j1 u1 … jq uq"，按 t 排序
        ll Tmax = atoll(argv[4]); ifstream f(argv[3]); string tag; vector<int> Nk(I.nk), Mj(I.nt);
        f >> tag; for (auto& x : Nk) f >> x; f >> tag; for (auto& x : Mj) f >> x;
        struct A { ll t; int s, r, rs; vector<pair<int,int>> tl; bool used = false; };
        vector<A> as; { ll t; int s2, r, rs, q; while (f >> t >> s2 >> r >> rs >> q) { A a{t, s2, r, rs, {}}; for (int c = 0; c < q; c++) { int j, u; f >> j >> u; a.tl.push_back({j, u}); } as.push_back(a); } }
        vector<int> rk, rph; vector<ll> rv;
        for (int k = 0; k < I.nk; k++) for (int c = 0; c < Nk[k]; c++) { rk.push_back(k); rph.push_back(0); rv.push_back(I.T[k].v0); }
        vector<vector<ll>> tf(I.nt); for (int j = 0; j < I.nt; j++) tf[j].assign(Mj[j], -1000000);
        int c0 = 0; for (int k = 0; k < I.nk; k++) c0 = max(c0, I.T[k].v0);
        const int np = I.psite.size(); vector<ll> fr(np); for (int q = 0; q < np; q++) fr[q] = c0 + I.S[I.psite[q]].lead;
        ll gaps = 0, duties = 0, resets = 0; map<string, ll> viol; size_t ptr = 0;
        ll t = *min_element(fr.begin(), fr.end());
        while (t < Tmax) {
            while (ptr < as.size() && as[ptr].t < t) { if (!as[ptr].used) viol["unmatched"]++; ptr++; }   // 没有岗位到点时的多余值岗
            vector<int> due; for (int q = 0; q < np; q++) if (fr[q] <= t) due.push_back(q);
            stable_sort(due.begin(), due.end(), [&](int a, int b) { return I.nelig[I.psite[a]] < I.nelig[I.psite[b]]; });
            for (int pq : due) {
                int s = I.psite[pq]; A* a = nullptr;
                for (size_t x = ptr; x < as.size() && as[x].t == t; x++) if (!as[x].used && as[x].s == s) { a = &as[x]; break; }
                if (!a) { gaps++; fr[pq] = t + 1; continue; }
                a->used = true; int r = a->r;
                if (r < 0 || r >= (int)rk.size()) { viol["robot_id"]++; fr[pq] = t + 1; gaps++; continue; }
                int k = rk[r], e = I.eidx[k][s]; const Type& T = I.T[k];
                if (e < 0) { viol["not_allowed"]++; fr[pq] = t + 1; gaps++; continue; }
                if (rv[r] + I.S[s].lead > t) viol["robot_not_ready"]++;
                int i = rph[r];
                if (a->rs && i == T.L - 1) viol["reset_at_last_phase"]++;
                ll need = t - T.pre[e]; double unin = 0; set<pair<int,int>> seenu;
                if ((int)a->tl.size() != T.q) viol["tool_count"]++;
                for (auto& ju : a->tl) {
                    int j = ju.first, u = ju.second;
                    if (find(T.C.begin(), T.C.end(), j) == T.C.end() || u < 0 || u >= Mj[j]) { viol["tool_type"]++; continue; }
                    if (!seenu.insert(ju).second) viol["tool_dup"]++;
                    if (tf[j][u] > need) viol["tool_busy"]++;
                    unin += I.TT[j].unin;
                }
                for (auto& ju : a->tl) { int j = ju.first, u = ju.second; if (j >= 0 && j < I.nt && u >= 0 && u < Mj[j]) tf[j][u] = t + T.d + postu(I, k, s, j, unin); }
                int ld = T.step[i][e].first, nx = T.step[i][e].second, nph = (i + 1) % T.L;
                if (a->rs) { nx = T.step[T.L - 1][e].second; nph = 0; resets++; }
                rv[r] = t - ld + nx; rph[r] = nph; duties++;
                fr[pq] = t + T.d;
            }
            for (size_t x = ptr; x < as.size() && as[x].t == t; x++) if (!as[x].used) { viol["unmatched"]++; as[x].used = true; }
            t = *min_element(fr.begin(), fr.end());
        }
        ll nv = 0; string vs = "{"; for (auto& kv : viol) { nv += kv.second; vs += (vs.size() > 1 ? ", \"" : "\"") + kv.first + "\": " + to_string(kv.second); } vs += "}";
        printf("{\"gaps\": %lld, \"violations\": %lld, \"viol\": %s, \"duties\": %lld, \"resets\": %lld, \"t_end\": %lld, \"n_assign\": %zu}\n", gaps, nv, vs.c_str(), duties, resets, t, as.size());
        return 0;
    }
    if (cmd == "rlp") {   // R_LP：定台数 + 复核
        RLP RR = parse_rlp(I, argv[3]); ll Hs = atoll(argv[4]), Tcert = atoll(argv[5]), Tver = atoll(argv[6]);
        Policy P; P.w.assign(13, 0); P.rlp = &RR; P.curdef = true;
        for (int a = 7; a < argc; a++) { string s = argv[a]; if (s.rfind("planf=", 0) == 0) { P.planf = parse_planf(I, s.substr(6)); P.planN = g_planN; } else if (s == "certadd=why") g_certadd = 0; else if (s == "certadd=trial60") g_certadd = 2;
            else if (s.rfind("sizelog=", 0) == 0) g_sizelog = fopen(s.substr(8).c_str(), "w");
            else if (s.rfind("planW=", 0) == 0) P.plan_until = atoll(s.c_str() + 6); }
        auto t0 = chrono::steady_clock::now();
        // 定台数：从 LP 取整出发跑 rost4 的"补足 / 回删 + 证书"循环；失败（加 60 次仍有缺岗）时把起点放大到 1.15、1.3、1.5 倍再跑，取第一个成功的
        const int facs[4] = {100, 115, 130, 150}; int fac_used = -1; SizeRes sr;
        for (int fi = 0; fi < 4; fi++) {
            g_N0.assign(I.nk, 0); g_M0.assign(I.nt, 0);
            for (int k = 0; k < I.nk; k++) g_N0[k] = (RR.N0[k] * facs[fi] + 99) / 100;
            for (int j = 0; j < I.nt; j++) g_M0[j] = (RR.M0[j] * facs[fi] + 99) / 100;
            sr = size_res(I, P, Hs, "certs", Tcert);
            if (sr.R >= 0) { fac_used = facs[fi]; break; }
        }
        double sec_size = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        Res v = (sr.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, sr.Nk, sr.Mj, P, Tver, true);
        int exact = (v.cert == 1) ? (exact_check(I, sr.Nk, sr.Mj, P, v.enter, v.period) ? 1 : 0) : -1;
        double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        int NR = 0, NT = 0; for (int x : sr.Nk) NR += x; for (int x : sr.Mj) NT += x;
        printf("{\"inst\": \"%s\", \"R\": %d, \"robots\": %d, \"tools\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, \"N0\": %s, \"M0\": %s, "
               "\"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"verify_enter\": %lld, \"verify_t_end\": %lld, \"exact\": %d, \"resets\": %lld, \"duties\": %lld, \"cyc_resets\": %lld, \"last_reset\": %lld, \"start_factor\": %d, \"sec_size\": %.2f, \"sec\": %.2f}\n",
               I.name.c_str(), sr.R, NR, NT, sr.status.c_str(), vec_s(sr.Nk).c_str(), vec_s(sr.Mj).c_str(), vec_s(RR.N0).c_str(), vec_s(RR.M0).c_str(),
               v.gaps, v.cert, v.period, v.enter, v.t_end, exact, v.resets, v.duties, v.cyc_resets, v.last_reset, fac_used, sec_size, sec);
        return 0;
    }
    if (cmd == "size") {   // T7：给定规则（可带 rdd=基因型文件）定台数 + Tverify 格复核
        ll Hs = atoll(argv[3]), Tcert = atoll(argv[4]), Tver = atoll(argv[5]);
        Rule U = parse_rule(I, argc, argv, 6); RDD DD;
        if (!U.geno.empty()) { Geno G = read_geno(U.geno); DD = decode_geno(G.g, G.W, G.ndev, G.qon, G.reset_opt); U.P.rdd = &DD; }
        auto t0 = chrono::steady_clock::now(); int fac = 0;
        SizeRes sr = size_rule(I, U.P, U.is_rlp, Hs, Tcert, &fac, U.facloop);
        Res v = (sr.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, sr.Nk, sr.Mj, U.P, Tver, true);
        int exact = (v.cert == 1) ? (exact_check(I, sr.Nk, sr.Mj, U.P, v.enter, v.period) ? 1 : 0) : -1;
        ll nm = sr.R < 0 ? -1 : near_miss(I, U.P, sr, Hs);
        double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        printf("{\"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, \"start_factor\": %d, \"near_miss\": %lld, \"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, \"last_reset\": %lld, \"ndec_w\": %lld, \"sec\": %.2f}\n",
               sr.R, sr.status.c_str(), vec_s(sr.Nk).c_str(), vec_s(sr.Mj).c_str(), fac, nm, v.gaps, v.cert, v.period, exact, v.resets, v.last_reset, v.ndec_w, sec);
        return 0;
    }
    if (cmd == "judge") {   // T11：固定台数完整判定。rost5 judge <inst> <Tcert> <Tverify> <Nk 逗号分隔> <Mj 逗号分隔> <规则选项…> [rdd=基因型文件]
        ll Tcert = atoll(argv[3]), Tver = atoll(argv[4]);
        auto csv = [](const string& x) { vector<int> v; stringstream ss(x); string y; while (getline(ss, y, ',')) if (!y.empty()) v.push_back(atoi(y.c_str())); return v; };
        vector<int> Nk = csv(argv[5]), Mj = csv(argv[6]);
        Rule U = parse_rule(I, argc, argv, 7); RDD DD;
        if (!U.geno.empty()) { Geno G = read_geno(U.geno); DD = decode_geno(G.g, G.W, G.ndev, G.qon, G.reset_opt); U.P.rdd = &DD; }
        auto t0 = chrono::steady_clock::now();
        Judge J = judge_fleet(I, Nk, Mj, U.P, Tcert, Tver);
        int R = 0; for (int x : Nk) R += x; for (int x : Mj) R += x;
        printf("{\"R\": %d, \"feasible\": %d, \"status\": \"%s\", \"gaps\": %lld, \"cert\": %d, \"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, \"last_reset\": %lld, \"sec\": %.2f}\n",
               R, J.feas ? 1 : 0, J.status.c_str(), J.gaps, J.cert, J.vgaps, J.vcert, J.vperiod, J.exact, J.resets, J.last_reset, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        return 0;
    }
    if (cmd == "fd") {   // T13：通用的定台数下降评价（与 rddfd 相同的框架，表示可换）。见 报告/T13_对手实现.md
        // rost5 fd <inst> <Hs> <Tcert> <Tverify> <Teval> <代数> <限时秒> <μ> <λ> <种子> <表示 rdd|w|rk|d5> <ρ> <q_on> <规则选项…> W=<格> [stall= swap= fdsel= fddump=]
        // 表示：rdd = RDD 的槽（与 rddfd 逐位相同）；w = 派工权重（ES：初始 3 个基础个体 + 正态随机，每代所有权重加 0.5·N(0,1)·(|x|+0.2)）；
        //       rk / d5 = 窗口内每个决策一个基因（初始个体 0 全 0，其余均匀随机；每个基因以 1/N_dev 的概率重抽）。
        ll Hs = atoll(argv[3]), Tcert = atoll(argv[4]), Tver = atoll(argv[5]), Teval = atoll(argv[6]);
        int gens = atoi(argv[7]); double tlim = atof(argv[8]); int mu = atoi(argv[9]), lam = atoi(argv[10]); ull seed = strtoull(argv[11], 0, 10);
        string repn = argv[12]; double rho = atof(argv[13]), qon = atof(argv[14]);
        Rule U = parse_rule(I, argc, argv, 15);
        const int RD = repn == "rdd" ? 0 : repn == "w" ? 1 : repn == "rk" ? 2 : repn == "d5" ? 3 : -1;
        if (RD < 0) { fprintf(stderr, "bad rep %s\n", repn.c_str()); return 1; }
        if (RD == 2 && U.P.pd_mode != 1) { U.P.pd_mode = 1; U.P.curdef = true; U.P.pdW = U.W; }
        if (RD == 3 && U.P.pd_mode != 2) { U.P.pd_mode = 2; U.P.repair = true; U.P.curdef = true; U.P.pdW = U.W; }
        vector<int> N0s = g_N0, M0s = g_M0;
        auto t0 = chrono::steady_clock::now();
        auto el = [&]() { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
        mt19937_64 rng(seed); uniform_real_distribution<double> U01(0, 1); normal_distribution<double> N01(0, 1);
        const int NW = (int)U.P.w.size();
        // 基因型 → 策略。Dec 持有 RDD 与基因的存储，不能复制
        struct Dec { Policy P; RDD D; vector<double> g; };
        int ndev = 0, K = 1;
        auto build = [&](Dec& d, const vector<double>& g) {
            d.P = U.P; d.g = g;
            if (RD == 0) { d.D = decode_geno(d.g, U.W, ndev, qon, U.reset_opt); d.P.rdd = &d.D; }
            else if (RD == 1) { d.P.w = d.g; }
            else { d.P.pdg = &d.g; }
        };
        vector<vector<double>> wbase = {{1, 0, 0, 0, 0, 0, 1, 0}, {0.1, 0, 10, 1, 0, 0, 1, 0.1}, {1, 0.5, 1, 0.5, 0, 0, 1, -0.1}};
        auto base_w = [&](int a) { vector<double> w(NW, 0.0); for (int q = 0; q < 8; q++) w[q] = wbase[a][q]; if (!U.P.lpz.empty()) { w[8] = -1; w[13] = 1.0; } return w; };
        // 1) F0：rdd / rk / d5 用默认规则（全 0 基因）定台数；w 用 3 个基础个体分别定台数，取资源最少的可行者（同资源取靠前的）
        int fac0 = 0; SizeRes s0; Res v0 = Res{-1, -2, -1, -1, {}, {}}; int ex0 = -1; bool base_feas = false; int base_hsfb = 0; int w_base_used = -1;
        vector<double> g0;   // F0 对应的基因型
        if (RD == 1) {
            for (int a = 0; a < 3; a++) {
                Dec d; build(d, base_w(a)); g_N0 = N0s; g_M0 = M0s; int fa = 0;
                SizeRes sa = size_rule(I, d.P, U.is_rlp, Hs, Tcert, &fa, U.facloop); g_N0 = N0s; g_M0 = M0s;
                if (a == 0) { s0 = sa; fac0 = fa; }
                if (sa.R < 0) continue;
                Res va = run(I, sa.Nk, sa.Mj, d.P, Tver, true); int ea = (va.cert == 1) ? (exact_check(I, sa.Nk, sa.Mj, d.P, va.enter, va.period) ? 1 : 0) : -1;
                if (!((va.cert == 1 && ea == 1) || va.gaps == 0)) continue;
                if (!base_feas || sa.R < s0.R) { s0 = sa; fac0 = fa; v0 = va; ex0 = ea; base_feas = true; w_base_used = a; g0 = d.g; }
            }
        } else {
            Dec d; build(d, {}); s0 = size_rule(I, d.P, U.is_rlp, Hs, Tcert, &fac0, U.facloop); base_hsfb = g_hsfb_used; g_N0 = N0s; g_M0 = M0s;
            if (s0.R >= 0) { v0 = run(I, s0.Nk, s0.Mj, d.P, Tver, true); ex0 = (v0.cert == 1) ? (exact_check(I, s0.Nk, s0.Mj, d.P, v0.enter, v0.period) ? 1 : 0) : -1; }
            base_feas = s0.R >= 0 && ((v0.cert == 1 && ex0 == 1) || v0.gaps == 0);
            // N_dev：默认规则在 F0（定不出台数时用最后试的台数）下窗口内的决策数
            if (!s0.Nk.empty()) {
                if (RD == 0) { RDD E; E.W = U.W; E.ndev = 0; E.maxpos = -1; Policy Pc = U.P; Pc.rdd = &E; ndev = (int)run(I, s0.Nk, s0.Mj, Pc, U.W, false).ndec_w; }
                else ndev = (int)run(I, s0.Nk, s0.Mj, d.P, U.W, false).ndec_w;
            }
            K = max(1, (int)ceil(rho * (double)ndev / qon - 1e-9));
        }
        const int glen = RD == 0 ? 2 * K : RD == 1 ? NW : max(ndev, 1);
        if (RD != 1) g0.assign(glen, 0.0);
        double sec0 = el();
        const int nres = I.nk + I.nt;
        vector<int> F(nres), lastF(nres);
        for (int k = 0; k < I.nk; k++) F[k] = s0.Nk.empty() ? 0 : s0.Nk[k];
        for (int j = 0; j < I.nt; j++) F[I.nk + j] = s0.Mj.empty() ? 0 : s0.Mj[j];
        lastF = F;
        auto split = [&](const vector<int>& f, vector<int>& Nk, vector<int>& Mj) { Nk.assign(f.begin(), f.begin() + I.nk); Mj.assign(f.begin() + I.nk, f.end()); };
        struct Ind { vector<double> g; vector<int> fl; ll f = 0; bool judged = false; };
        long long nevals = 0, nreuse = 0, njudge = 0, njudge_fail = 0, nquick = 0, ndec_steps = 0, nback = 0; double sec_eval = 0, sec_judge = 0, sec_dec = 0;
        FILE* dump = g_fddump.empty() ? nullptr : fopen(g_fddump.c_str(), "w"); int ndump = 0;
        auto gaps_at = [&](const vector<double>& g, const vector<int>& f, ll cap) -> ll {
            Dec d; build(d, g); vector<int> Nk, Mj; split(f, Nk, Mj);
            ll r = run(I, Nk, Mj, d.P, Teval, false, -1, nullptr, cap).gaps;
            if (dump && ndump < 40 && (nevals % 37 == 5)) {
                ll ex = run(I, Nk, Mj, d.P, Teval, false).gaps; ndump++;
                fprintf(dump, "%s\t%s\t%lld\t", vec_s(Nk).c_str(), vec_s(Mj).c_str(), ex);
                for (size_t a = 0; a < g.size(); a++) fprintf(dump, "%s%.17g", a ? "," : "", g[a]);
                fprintf(dump, "\n"); fflush(dump);
            }
            return r;
        };
        auto evalx = [&](Ind& x, ll cap) { auto a = chrono::steady_clock::now(); x.f = gaps_at(x.g, x.fl, cap); x.judged = false; nevals++;
                                           sec_eval += chrono::duration<double>(chrono::steady_clock::now() - a).count(); };
        auto choose_dec = [&](const vector<double>& g, const set<int>& tabu) -> int {
            auto a = chrono::steady_clock::now(); ll bg = LLONG_MAX; int bi = -1; vector<int> f = lastF;
            for (int idx = 0; idx < nres; idx++) {
                if (f[idx] == 0 || tabu.count(idx)) continue;
                f[idx]--; ll gg = gaps_at(g, f, bg); f[idx]++;
                if (bi < 0 || gg < bg) { bg = gg; bi = idx; }
            }
            sec_dec += chrono::duration<double>(chrono::steady_clock::now() - a).count(); ndec_steps++;
            return bi;
        };
        // 2) 初始种群：rdd 与 rddfd 相同（个体 0 全关，其余均匀随机）；rk / d5 个体 0 全 0，其余均匀随机；w 个体 0–2 为基础个体，其余正态随机（与 es 命令相同）
        vector<Ind> pop;
        for (int a = 0; a < mu + lam; a++) {
            Ind x;
            if (RD == 1) { if (a < 3) x.g = base_w(a); else { x.g.assign(NW, 0.0); for (auto& v : x.g) v = N01(rng); } }
            else { x.g.assign(glen, 0.0); if (a > 0) for (auto& v : x.g) v = U01(rng); }
            pop.push_back(x);
        }
        int bestR = base_feas ? s0.R : -1; vector<double> bestg = g0; vector<int> bestF = F; string best_status = s0.status;
        Judge bestJ; bestJ.feas = base_feas; bestJ.status = s0.status; bestJ.vgaps = v0.gaps; bestJ.vcert = v0.cert; bestJ.vperiod = v0.period; bestJ.exact = ex0; bestJ.resets = v0.resets; bestJ.last_reset = v0.last_reset;
        int f0_src = base_feas ? 0 : -1; double sec_f0 = 0;
        if (!base_feas) {   // 默认（或基础个体）定不出台数：按顺序对初始种群的其余个体完整定台数，取第一个可行的
            auto a0 = chrono::steady_clock::now();
            for (int a = (RD == 1 ? 3 : 1); a < mu + lam && el() < tlim; a++) {
                Dec d; build(d, pop[a].g);
                g_N0 = N0s; g_M0 = M0s; SizeRes sa = size_rule(I, d.P, U.is_rlp, Hs, Tcert, nullptr, U.facloop); g_N0 = N0s; g_M0 = M0s;
                if (sa.R < 0) continue;
                Judge J = judge_fleet(I, sa.Nk, sa.Mj, d.P, Tcert, Tver);
                if (!J.feas) continue;
                for (int k = 0; k < I.nk; k++) F[k] = sa.Nk[k];
                for (int j = 0; j < I.nt; j++) F[I.nk + j] = sa.Mj[j];
                lastF = F; bestF = F; bestR = sa.R; bestg = pop[a].g; best_status = J.status; bestJ = J; base_feas = true; f0_src = a; break;
            }
            sec_f0 = chrono::duration<double>(chrono::steady_clock::now() - a0).count();
        }
        string traj = "["; char buf[256];
        auto add_traj = [&](int gi) { snprintf(buf, 256, "%s[%d, %d, %lld, %.1f, \"%s\"]", traj.size() > 1 ? ", " : "", bestR, gi, nevals, el(), best_status.c_str()); traj += buf; };
        add_traj(-1);
        set<int> tabu; int cur_r = -1; int stall = 0; int gi = 0; long long gens_done = 0;
        auto by_f = [](const Ind& a, const Ind& b) { return a.f < b.f; };
        auto set_target = [&](const vector<double>& g) -> bool {
            int r = choose_dec(g, tabu);
            if (r < 0 && !tabu.empty()) { tabu.clear(); r = choose_dec(g, tabu); }
            if (r < 0) return false;
            cur_r = r; F = lastF; F[r]--; stall = 0;
            for (auto& x : pop) { x.fl = F; evalx(x, LLONG_MAX); }
            stable_sort(pop.begin(), pop.end(), by_f);
            return true;
        };
        bool alive = base_feas;
        if (alive) { pop.resize(mu + lam); alive = set_target(g_fdsel == "base" ? g0 : bestg); pop.resize(mu); }
        const double pm = RD == 0 ? 1.0 / (2.0 * K) : 1.0 / (double)glen;
        while (alive && gi < gens && el() < tlim) {
            vector<Ind> kids; ll cap = pop.back().f;
            for (int c = 0; c < lam; c++) {
                const Ind& pa = pop[uniform_int_distribution<int>(0, mu - 1)(rng)];
                Ind x; x.g = pa.g; x.fl = pa.fl; bool ch = false;
                if (RD == 1) { for (auto& v : x.g) v = v + 0.5 * N01(rng) * (fabs(v) + 0.2); ch = true; }
                else for (auto& v : x.g) if (U01(rng) < pm) { v = U01(rng); ch = true; }
                if (U01(rng) < g_swap) {
                    vector<int> pos; for (int idx = 0; idx < nres; idx++) if (x.fl[idx] > 0) pos.push_back(idx);
                    if (!pos.empty() && nres > 1) {
                        int a = pos[uniform_int_distribution<int>(0, (int)pos.size() - 1)(rng)];
                        int b = uniform_int_distribution<int>(0, nres - 2)(rng); if (b >= a) b++;
                        x.fl[a]--; x.fl[b]++; ch = true;
                    }
                }
                if (!ch) { x.f = pa.f; x.judged = pa.judged; nreuse++; }
                else evalx(x, cap);
                kids.push_back(x);
            }
            for (auto& x : pop) kids.push_back(x);
            pop.swap(kids); stable_sort(pop.begin(), pop.end(), by_f); pop.resize(mu);
            gi++; gens_done++;
            bool found = false;
            for (size_t q = 0; q < pop.size() && !found; q++) {
                Ind& x = pop[q];
                if (x.f != 0 || x.judged) continue;
                auto a = chrono::steady_clock::now();
                Dec d; build(d, x.g); vector<int> Nk, Mj; split(x.fl, Nk, Mj);
                ll capJ = max(pop.back().f, (ll)0);
                ll qg = run(I, Nk, Mj, d.P, Tcert, false, -1, nullptr, capJ).gaps;
                if (qg > 0) { njudge++; njudge_fail++; nquick++; x.judged = true; x.f = qg;
                              sec_judge += chrono::duration<double>(chrono::steady_clock::now() - a).count(); continue; }
                Judge J = judge_fleet(I, Nk, Mj, d.P, Tcert, Tver); njudge++; x.judged = true;
                sec_judge += chrono::duration<double>(chrono::steady_clock::now() - a).count();
                if (J.feas) {
                    found = true; bestR = 0; for (int v : x.fl) bestR += v; bestg = x.g; bestF = x.fl; best_status = J.status; bestJ = J;
                    add_traj(gi); lastF = x.fl; tabu.clear();
                    alive = set_target(g_fdsel == "base" ? g0 : bestg);
                } else { njudge_fail++; x.f = J.gaps; }
            }
            if (!found) {
                stable_sort(pop.begin(), pop.end(), by_f);
                if (++stall >= g_stall) { tabu.insert(cur_r); nback++; alive = set_target(g_fdsel == "base" ? g0 : bestg); }
            }
        }
        if (dump) fclose(dump);
        int nact = 0;
        if (RD == 0) { RDD Db = decode_geno(bestg, U.W, ndev, qon, U.reset_opt); for (double u : Db.dev) if (u >= 0) nact++; }
        else if (RD >= 2) for (double u : bestg) if (u > 0) nact++;
        int nrep = -1; if (RD == 3 && bestR >= 0) { Dec d; build(d, bestg); vector<int> Nk, Mj; split(bestF, Nk, Mj); vector<int> f = d5_repair(I, Nk, Mj, d.P); nrep = 0; for (int v : f) if (v >= 0) nrep++; }
        traj += "]";
        string gs = "["; for (size_t i = 0; i < bestg.size(); i++) { snprintf(buf, 40, "%s%.17g", i ? "," : "", bestg[i]); gs += buf; } gs += "]";
        vector<int> bN, bM; split(bestF, bN, bM);
        printf("{\"rep\": \"%s\", \"base_R\": %d, \"base_status\": \"%s\", \"base_Nk\": %s, \"base_Mj\": %s, \"base_factor\": %d, \"base_hsfb\": %d, \"w_base_used\": %d, \"base_verify_gaps\": %lld, \"base_verify_cert\": %d, \"base_exact\": %d, \"base_sec\": %.2f, "
               "\"W\": %lld, \"ndev\": %d, \"K\": %d, \"glen\": %d, \"rho\": %g, \"qon\": %g, \"seed\": %llu, \"Teval\": %lld, \"fdsel\": \"%s\", \"stall\": %d, \"swap\": %g, "
               "\"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, \"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, \"last_reset\": %lld, \"base_resets\": %lld, \"base_last_reset\": %lld, "
               "\"f0_src\": %d, \"sec_f0\": %.2f, \"n_active\": %d, \"n_repair\": %d, \"gens\": %lld, \"evals\": %lld, \"reused\": %lld, \"judges\": %lld, \"judge_fail\": %lld, \"judge_quick\": %lld, \"dec_steps\": %lld, \"backtracks\": %lld, "
               "\"sec_eval\": %.2f, \"sec_judge\": %.2f, \"sec_dec\": %.2f, \"sec\": %.2f, \"traj\": %s, \"genes\": %s}\n",
               repn.c_str(), s0.R, s0.status.c_str(), vec_s(s0.Nk).c_str(), vec_s(s0.Mj).c_str(), fac0, base_hsfb, w_base_used, v0.gaps, v0.cert, ex0, sec0,
               U.W, ndev, K, glen, rho, qon, seed, Teval, g_fdsel.c_str(), g_stall, g_swap,
               bestR, best_status.c_str(), vec_s(bN).c_str(), vec_s(bM).c_str(), bestJ.vgaps, bestJ.vcert, bestJ.vperiod, bestJ.exact, bestJ.resets, bestJ.last_reset, v0.resets, v0.last_reset,
               f0_src, sec_f0, nact, nrep, gens_done, nevals, nreuse, njudge, njudge_fail, nquick, ndec_steps, nback, sec_eval, sec_judge, sec_dec, el(), traj.c_str(), gs.c_str());
        return 0;
    }
    if (cmd == "rddfd") {   // T11：RDD + 定台数下降评价（见 报告/T11_便宜评价.md）
        // rost5 rddfd <inst> <Hs> <Tcert> <Tverify> <Teval> <代数> <限时秒> <μ> <λ> <种子> <ρ> <q_on> <规则选项…> W=<格> [fdsel=best|base] [stall=代] [fddump=文件]
        ll Hs = atoll(argv[3]), Tcert = atoll(argv[4]), Tver = atoll(argv[5]), Teval = atoll(argv[6]);
        int gens = atoi(argv[7]); double tlim = atof(argv[8]); int mu = atoi(argv[9]), lam = atoi(argv[10]); ull seed = strtoull(argv[11], 0, 10);
        double rho = atof(argv[12]), qon = atof(argv[13]);
        Rule U = parse_rule(I, argc, argv, 14);
        vector<int> N0s = g_N0, M0s = g_M0;
        auto t0 = chrono::steady_clock::now();
        auto el = [&]() { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
        // 1) 默认规则本身（与 rdd 命令相同）：F0 = 默认规则定出的台数；N_dev、K 也相同
        int fac0 = 0; SizeRes s0 = size_rule(I, U.P, U.is_rlp, Hs, Tcert, &fac0, U.facloop); int base_hsfb = g_hsfb_used;
        g_N0 = N0s; g_M0 = M0s;
        Res v0 = (s0.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, s0.Nk, s0.Mj, U.P, Tver, true);
        int ex0 = (v0.cert == 1) ? (exact_check(I, s0.Nk, s0.Mj, U.P, v0.enter, v0.period) ? 1 : 0) : -1;
        bool base_feas = s0.R >= 0 && ((v0.cert == 1 && ex0 == 1) || v0.gaps == 0);
        RDD E; E.W = U.W; E.ndev = 0; E.maxpos = -1; Policy Pc = U.P; Pc.rdd = &E;
        int ndev = s0.Nk.empty() ? 0 : (int)run(I, s0.Nk, s0.Mj, Pc, U.W, false).ndec_w;   // 与 rdd 命令相同（默认规则定不出台数时用它最后试的台数）
        int K = max(1, (int)ceil(rho * (double)ndev / qon - 1e-9));
        double sec0 = el();
        const int nres = I.nk + I.nt;
        // 当前目标台数 F 与最近一个判定可行的解
        vector<int> F(nres), lastF(nres);
        for (int k = 0; k < I.nk; k++) F[k] = s0.Nk.empty() ? 0 : s0.Nk[k];
        for (int j = 0; j < I.nt; j++) F[I.nk + j] = s0.Mj.empty() ? 0 : s0.Mj[j];
        lastF = F;
        auto split = [&](const vector<int>& f, vector<int>& Nk, vector<int>& Mj) { Nk.assign(f.begin(), f.begin() + I.nk); Mj.assign(f.begin() + I.nk, f.end()); };
        struct Ind { vector<double> g; vector<int> fl; ll f = 0; bool judged = false; };   // fl：该个体的台数（总数 = 当前目标）
        long long nevals = 0, nreuse = 0, njudge = 0, njudge_fail = 0, nquick = 0, ndec_steps = 0, nback = 0; double sec_eval = 0, sec_judge = 0, sec_dec = 0;
        FILE* dump = g_fddump.empty() ? nullptr : fopen(g_fddump.c_str(), "w"); int ndump = 0;
        auto gaps_at = [&](const vector<double>& g, const vector<int>& f, ll cap) -> ll {
            RDD D = decode_geno(g, U.W, ndev, qon, U.reset_opt); Policy P = U.P; P.rdd = &D; vector<int> Nk, Mj; split(f, Nk, Mj);
            ll r = run(I, Nk, Mj, P, Teval, false, -1, nullptr, cap).gaps;
            if (dump && ndump < 40 && (nevals % 37 == 5)) {   // 抽样：精确缺岗（不设上限）写进核对文件
                ll ex = run(I, Nk, Mj, P, Teval, false).gaps; ndump++;
                fprintf(dump, "%s\t%s\t%lld\t", vec_s(Nk).c_str(), vec_s(Mj).c_str(), ex);
                for (size_t a = 0; a < g.size(); a++) fprintf(dump, "%s%.17g", a ? "," : "", g[a]);
                fprintf(dump, "\n"); fflush(dump);
            }
            return r;
        };
        auto evalx = [&](Ind& x, ll cap) { auto a = chrono::steady_clock::now(); x.f = gaps_at(x.g, x.fl, cap); x.judged = false; nevals++;
                                           sec_eval += chrono::duration<double>(chrono::steady_clock::now() - a).count(); };
        // 选减哪个资源：对每种资源 r 试 lastF − e_r，用基因型 g 跑 T_eval 格，取缺岗最少的；并列取下标小的；跳过已试过的（tabu）
        auto choose_dec = [&](const vector<double>& g, const set<int>& tabu) -> int {
            auto a = chrono::steady_clock::now(); ll bg = LLONG_MAX; int bi = -1; vector<int> f = lastF;
            for (int idx = 0; idx < nres; idx++) {
                if (f[idx] == 0 || tabu.count(idx)) continue;
                f[idx]--; ll gg = gaps_at(g, f, bg); f[idx]++;
                if (bi < 0 || gg < bg) { bg = gg; bi = idx; }
            }
            sec_dec += chrono::duration<double>(chrono::steady_clock::now() - a).count(); ndec_steps++;
            return bi;
        };
        // 2) 初始种群（与 rdd 命令相同的随机数流：个体 0 全关，其余均匀随机）
        mt19937_64 rng(seed); uniform_real_distribution<double> U01(0, 1);
        vector<Ind> pop;
        for (int a = 0; a < mu + lam; a++) { Ind x; x.g.assign(2 * K, 0.0); if (a > 0) for (auto& v : x.g) v = U01(rng); pop.push_back(x); }
        // 最好解：先是默认规则本身
        int bestR = base_feas ? s0.R : -1; vector<double> bestg(2 * K, 0.0); vector<int> bestF = F; string best_status = s0.status;
        Judge bestJ; bestJ.feas = base_feas; bestJ.status = s0.status; bestJ.vgaps = v0.gaps; bestJ.vcert = v0.cert; bestJ.vperiod = v0.period; bestJ.exact = ex0; bestJ.resets = v0.resets; bestJ.last_reset = v0.last_reset;
        // F0 的来源：0 = 默认规则；默认规则定不出台数时，按顺序对初始种群的个体 1, 2, … 完整定台数（与 rdd 命令的初始评价相同），
        // 取第一个定出台数且完整判定可行的个体 a（f0_src = a）；都不行时为 −1，搜索不进行。
        int f0_src = base_feas ? 0 : -1; double sec_f0 = 0;
        if (!base_feas) {
            auto a0 = chrono::steady_clock::now();
            for (int a = 1; a < mu + lam && el() < tlim; a++) {
                RDD D = decode_geno(pop[a].g, U.W, ndev, qon, U.reset_opt); Policy P = U.P; P.rdd = &D;
                g_N0 = N0s; g_M0 = M0s; SizeRes sa = size_rule(I, P, U.is_rlp, Hs, Tcert, nullptr, U.facloop); g_N0 = N0s; g_M0 = M0s;
                if (sa.R < 0) continue;
                Judge J = judge_fleet(I, sa.Nk, sa.Mj, P, Tcert, Tver);
                if (!J.feas) continue;
                for (int k = 0; k < I.nk; k++) F[k] = sa.Nk[k];
                for (int j = 0; j < I.nt; j++) F[I.nk + j] = sa.Mj[j];
                lastF = F; bestF = F; bestR = sa.R; bestg = pop[a].g; best_status = J.status; bestJ = J; base_feas = true; f0_src = a; break;
            }
            sec_f0 = chrono::duration<double>(chrono::steady_clock::now() - a0).count();
        }
        string traj = "[";
        char buf[256];
        auto add_traj = [&](int gi) { snprintf(buf, 256, "%s[%d, %d, %lld, %.1f, \"%s\"]", traj.size() > 1 ? ", " : "", bestR, gi, nevals, el(), best_status.c_str()); traj += buf; };
        add_traj(-1);
        set<int> tabu; int cur_r = -1; int stall = 0; int gi = 0; long long gens_done = 0;
        auto by_f = [](const Ind& a, const Ind& b) { return a.f < b.f; };
        auto set_target = [&](const vector<double>& g) -> bool {    // 从 lastF 选一个资源减 1，作为新目标；全部试过就清空 tabu 再选
            int r = choose_dec(g, tabu);
            if (r < 0 && !tabu.empty()) { tabu.clear(); r = choose_dec(g, tabu); }
            if (r < 0) return false;
            cur_r = r; F = lastF; F[r]--; stall = 0;
            for (auto& x : pop) { x.fl = F; evalx(x, LLONG_MAX); }
            stable_sort(pop.begin(), pop.end(), by_f);
            return true;
        };
        bool alive = base_feas;
        if (alive) {
            pop.resize(mu + lam);
            alive = set_target(g_fdsel == "base" ? vector<double>(2 * K, 0.0) : bestg);
            pop.resize(mu);
        }
        const double pm = 1.0 / (2.0 * K);
        while (alive && gi < gens && el() < tlim) {
            // 子代：与 rdd 命令相同的变异
            vector<Ind> kids; ll cap = pop.back().f;
            for (int c = 0; c < lam; c++) {
                const Ind& pa = pop[uniform_int_distribution<int>(0, mu - 1)(rng)];
                Ind x; x.g = pa.g; x.fl = pa.fl; bool ch = false;
                for (auto& v : x.g) if (U01(rng) < pm) { v = U01(rng); ch = true; }
                if (U01(rng) < g_swap) {   // 台数互换：从台数 > 0 的资源里均匀取一个 −1，再从其余资源里均匀取一个 +1
                    vector<int> pos; for (int idx = 0; idx < nres; idx++) if (x.fl[idx] > 0) pos.push_back(idx);
                    if (!pos.empty() && nres > 1) {
                        int a = pos[uniform_int_distribution<int>(0, (int)pos.size() - 1)(rng)];
                        int b = uniform_int_distribution<int>(0, nres - 2)(rng); if (b >= a) b++;
                        x.fl[a]--; x.fl[b]++; ch = true;
                    }
                }
                if (!ch) { x.f = pa.f; x.judged = pa.judged; nreuse++; }
                else evalx(x, cap);
                kids.push_back(x);
            }
            // 选择：子代在前（同缺岗时子代胜出，允许中性漂移），稳定排序取前 μ
            for (auto& x : pop) kids.push_back(x);
            pop.swap(kids); stable_sort(pop.begin(), pop.end(), by_f); pop.resize(mu);
            gi++; gens_done++;
            // 判定：T_eval 内 0 缺岗且没判定过的个体，按顺序做完整判定
            bool found = false;
            for (size_t q = 0; q < pop.size() && !found; q++) {
                Ind& x = pop[q];
                if (x.f != 0 || x.judged) continue;
                auto a = chrono::steady_clock::now();
                RDD D = decode_geno(x.g, U.W, ndev, qon, U.reset_opt); Policy P = U.P; P.rdd = &D; vector<int> Nk, Mj; split(x.fl, Nk, Mj);
                // 快速否决：先不判环跑 Tcert 格，缺岗超过种群里最差的适应度就停。这里有缺岗 ⇔ 完整判定的判环阶段有缺岗（仿真是确定性的，
                // 判出无缺岗的循环就不会再有缺岗），所以否决与完整判定的结论相同；不可行个体的适应度取这次的缺岗数（超过上限时为上限 + 1）。
                ll capJ = max(pop.back().f, (ll)0);
                ll qg = run(I, Nk, Mj, P, Tcert, false, -1, nullptr, capJ).gaps;
                if (qg > 0) { njudge++; njudge_fail++; nquick++; x.judged = true; x.f = qg;
                              sec_judge += chrono::duration<double>(chrono::steady_clock::now() - a).count(); continue; }
                Judge J = judge_fleet(I, Nk, Mj, P, Tcert, Tver); njudge++; x.judged = true;
                sec_judge += chrono::duration<double>(chrono::steady_clock::now() - a).count();
                if (J.feas) {
                    found = true; bestR = 0; for (int v : x.fl) bestR += v; bestg = x.g; bestF = x.fl; best_status = J.status; bestJ = J;
                    add_traj(gi); lastF = x.fl; tabu.clear();
                    vector<double> gsel = g_fdsel == "base" ? vector<double>(2 * K, 0.0) : bestg;
                    alive = set_target(gsel);
                } else { njudge_fail++; x.f = J.gaps; }
            }
            if (!found) {
                stable_sort(pop.begin(), pop.end(), by_f);
                if (++stall >= g_stall) { tabu.insert(cur_r); nback++; alive = set_target(g_fdsel == "base" ? vector<double>(2 * K, 0.0) : bestg); }
            }
        }
        if (dump) fclose(dump);
        RDD Db = decode_geno(bestg, U.W, ndev, qon, U.reset_opt); int nact = 0; for (double u : Db.dev) if (u >= 0) nact++;
        traj += "]";
        string gs = "["; for (size_t i = 0; i < bestg.size(); i++) { snprintf(buf, 40, "%s%.17g", i ? "," : "", bestg[i]); gs += buf; } gs += "]";
        vector<int> bN, bM; split(bestF, bN, bM);
        printf("{\"base_R\": %d, \"base_status\": \"%s\", \"base_Nk\": %s, \"base_Mj\": %s, \"base_factor\": %d, \"base_hsfb\": %d, \"base_verify_gaps\": %lld, \"base_verify_cert\": %d, \"base_exact\": %d, \"base_sec\": %.2f, "
               "\"W\": %lld, \"ndev\": %d, \"K\": %d, \"rho\": %g, \"qon\": %g, \"seed\": %llu, \"Teval\": %lld, \"fdsel\": \"%s\", \"stall\": %d, \"swap\": %g, "
               "\"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, \"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, \"last_reset\": %lld, \"base_resets\": %lld, \"base_last_reset\": %lld, "
               "\"f0_src\": %d, \"sec_f0\": %.2f, \"n_active\": %d, \"gens\": %lld, \"evals\": %lld, \"reused\": %lld, \"judges\": %lld, \"judge_fail\": %lld, \"judge_quick\": %lld, \"dec_steps\": %lld, \"backtracks\": %lld, "
               "\"sec_eval\": %.2f, \"sec_judge\": %.2f, \"sec_dec\": %.2f, \"sec\": %.2f, \"traj\": %s, \"genes\": %s}\n",
               s0.R, s0.status.c_str(), vec_s(s0.Nk).c_str(), vec_s(s0.Mj).c_str(), fac0, base_hsfb, v0.gaps, v0.cert, ex0, sec0,
               U.W, ndev, K, rho, qon, seed, Teval, g_fdsel.c_str(), g_stall, g_swap,
               bestR, best_status.c_str(), vec_s(bN).c_str(), vec_s(bM).c_str(), bestJ.vgaps, bestJ.vcert, bestJ.vperiod, bestJ.exact, bestJ.resets, bestJ.last_reset, v0.resets, v0.last_reset,
               f0_src, sec_f0, nact, gens_done, nevals, nreuse, njudge, njudge_fail, nquick, ndec_steps, nback, sec_eval, sec_judge, sec_dec, el(), traj.c_str(), gs.c_str());
        return 0;
    }
    if (cmd == "rdd") {    // T7/T8：RDD 的 (μ + λ) 进化。基因型 K 个槽 (p, a)，每个基因以 1/(2K) 的概率重抽
        ll Hs = atoll(argv[3]), Tcert = atoll(argv[4]), Tver = atoll(argv[5]);
        int gens = atoi(argv[6]), mu = atoi(argv[7]), lam = atoi(argv[8]); ull seed = strtoull(argv[9], 0, 10);
        double rho = atof(argv[10]), qon = atof(argv[11]);
        Rule U = parse_rule(I, argc, argv, 12);
        vector<int> N0s = g_N0, M0s = g_M0;     // 线性规则的定台数起点（size_rule 对 R_LP 会改写 g_N0，这里线性规则不受影响）
        auto t0 = chrono::steady_clock::now();
        // 1) 默认规则本身：定台数 + 复核；N_dev = 默认规则在其台数下窗口 [0, W) 内的决策数
        int fac0 = 0; SizeRes s0 = size_rule(I, U.P, U.is_rlp, Hs, Tcert, &fac0, U.facloop);
        g_N0 = N0s; g_M0 = M0s;
        Res v0 = (s0.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, s0.Nk, s0.Mj, U.P, Tver, true);
        int ex0 = (v0.cert == 1) ? (exact_check(I, s0.Nk, s0.Mj, U.P, v0.enter, v0.period) ? 1 : 0) : -1;
        RDD E; E.W = U.W; E.ndev = 0; E.maxpos = -1; Policy Pc = U.P; Pc.rdd = &E;
        int ndev = (int)run(I, s0.Nk, s0.Mj, Pc, U.W, false).ndec_w;
        int K = max(1, (int)ceil(rho * (double)ndev / qon - 1e-9));
        double sec0 = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        // 2) 进化
        mt19937_64 rng(seed); uniform_real_distribution<double> U01(0, 1);
        struct Ind { vector<double> g; double f1; ll f2; SizeRes s; };
        long long ndec_run = 0;
        auto fit = [&](Ind& x) {
            RDD D = decode_geno(x.g, U.W, ndev, qon, U.reset_opt); Policy P = U.P; P.rdd = &D;
            g_N0 = N0s; g_M0 = M0s;
            x.s = size_rule(I, P, U.is_rlp, Hs, Tcert, nullptr, U.facloop); ndec_run++;
            x.f1 = x.s.R < 0 ? 999 : x.s.R + (x.s.status == "uncert" ? 3.0 : 0.0);
            x.f2 = x.s.R < 0 ? LLONG_MAX : near_miss(I, P, x.s, Hs);
        };
        auto better = [](const Ind& a, const Ind& b) { return a.f1 < b.f1 || (a.f1 == b.f1 && a.f2 < b.f2); };
        vector<Ind> pop;
        for (int a = 0; a < mu + lam; a++) {
            Ind x; x.g.assign(2 * K, 0.0);
            if (a > 0) for (auto& v : x.g) v = U01(rng);      // 个体 0：全部槽关闭（即默认规则）
            fit(x); pop.push_back(x);
        }
        stable_sort(pop.begin(), pop.end(), better); pop.resize(mu);
        vector<double> hist = {pop[0].f1}; long long nsame = 0;
        const double pm = 1.0 / (2.0 * K);
        for (int gi = 0; gi < gens; gi++) {
            vector<Ind> kids;
            for (int c = 0; c < lam; c++) {
                const Ind& pa = pop[uniform_int_distribution<int>(0, mu - 1)(rng)];
                Ind x; x.g = pa.g; bool ch = false;
                for (auto& v : x.g) if (U01(rng) < pm) { v = U01(rng); ch = true; }
                if (!ch) { x.f1 = pa.f1; x.f2 = pa.f2; x.s = pa.s; nsame++; }   // 没有基因变化：解码结果相同，直接沿用
                else fit(x);
                kids.push_back(x);
            }
            for (auto& k : kids) pop.push_back(k);
            stable_sort(pop.begin(), pop.end(), better); pop.resize(mu);
            hist.push_back(pop[0].f1);
        }
        const Ind& b = pop[0];
        RDD Db = decode_geno(b.g, U.W, ndev, qon, U.reset_opt); Policy Pb = U.P; Pb.rdd = &Db;
        Res v = (b.s.R < 0) ? Res{-1, -2, -1, -1, {}, {}} : run(I, b.s.Nk, b.s.Mj, Pb, Tver, true);
        int exact = (v.cert == 1) ? (exact_check(I, b.s.Nk, b.s.Mj, Pb, v.enter, v.period) ? 1 : 0) : -1;
        int nact = 0; for (double u : Db.dev) if (u >= 0) nact++;
        double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        string hs = "["; char buf[40]; for (size_t i = 0; i < hist.size(); i++) { snprintf(buf, 40, "%s%.1f", i ? "," : "", hist[i]); hs += buf; } hs += "]";
        string gs = "["; for (size_t i = 0; i < b.g.size(); i++) { snprintf(buf, 40, "%s%.17g", i ? "," : "", b.g[i]); gs += buf; } gs += "]";
        printf("{\"base_R\": %d, \"base_status\": \"%s\", \"base_Nk\": %s, \"base_Mj\": %s, \"base_factor\": %d, \"base_verify_gaps\": %lld, \"base_verify_cert\": %d, \"base_exact\": %d, \"base_sec\": %.2f, "
               "\"W\": %lld, \"ndev\": %d, \"K\": %d, \"rho\": %g, \"qon\": %g, \"seed\": %llu, "
               "\"R\": %d, \"status\": \"%s\", \"Nk\": %s, \"Mj\": %s, \"near_miss\": %lld, \"verify_gaps\": %lld, \"verify_cert\": %d, \"verify_period\": %lld, \"exact\": %d, \"resets\": %lld, "
               "\"n_active\": %d, \"decodes\": %lld, \"reused\": %lld, \"sec\": %.1f, \"hist\": %s, \"genes\": %s}\n",
               s0.R, s0.status.c_str(), vec_s(s0.Nk).c_str(), vec_s(s0.Mj).c_str(), fac0, v0.gaps, v0.cert, ex0, sec0,
               U.W, ndev, K, rho, qon, seed,
               b.s.R, b.s.status.c_str(), vec_s(b.s.Nk).c_str(), vec_s(b.s.Mj).c_str(), b.f2 == LLONG_MAX ? -1 : b.f2, v.gaps, v.cert, v.period, exact, v.resets,
               nact, ndec_run, nsame, sec, hs.c_str(), gs.c_str());
        return 0;
    }
    if (cmd == "es") {
        string mode = argv[3]; ll H = atoll(argv[4]); int gens = atoi(argv[5]), mu = atoi(argv[6]), lam = atoi(argv[7]);
        ull seed = strtoull(argv[8], 0, 10); ll Tcert = atoll(argv[9]), Tver = atoll(argv[10]);
        bool gp = false, o_reset = false, o_curdef = false, o_hard = false, o_pad = false, o_stag = false; ll o_Tw = 0; vector<vector<int>> o_quota; vector<vector<vector<pair<ll,int>>>> o_plan; vector<vector<vector<pair<ll,double>>>> o_planf; vector<vector<double>> o_lpz;
        for (int a = 11; a < argc; a++) { string s = argv[a];
            if (s == "gp") gp = true; else if (s == "reset") o_reset = true; else if (s == "curdef") o_curdef = true; else if (s == "hard") o_hard = true; else if (s == "pad") o_pad = true; else if (s == "stagger") o_stag = true; else if (s.rfind("quota=", 0) == 0) o_quota = parse_quota(I, s.substr(6)); else if (s.rfind("plan=", 0) == 0) o_plan = parse_plan(I, s.substr(5)); else if (s.rfind("planf=", 0) == 0) o_planf = parse_planf(I, s.substr(6)); else if (s.rfind("lpg=", 0) == 0) o_lpz = parse_lpz(I, s.substr(4), g_N0, g_M0);
            else if (s.rfind("Tw=", 0) == 0) o_Tw = atoll(s.c_str() + 3); else if (s == "certadd=why") g_certadd = 0; else if (s == "certadd=trial60") g_certadd = 2; }
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
