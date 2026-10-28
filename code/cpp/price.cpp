// 列生成定价：单台 k 型机器人在时间展开 DAG 上的最长路（收集的对偶收益最大）。
// 节点 (相位 i, 虚拟时刻 v)，v ∈ [v0, H)。弧：
//   等待 (i, v) → (i, v+1)；
//   值岗 (i, v) →(站 e) ((i+1) mod L, v + nx[i][e])，收益 w[e][v]；
//   提前大修（reset=1 且 i ≠ L−1）(i, v) →(站 e) (0, v + nxR[e])，收益同为 w[e][v]。
// 弧头 ≥ H 即进入汇点。起点 (0, v0)，值为 0。
// 输出：最大收益、路径上的值岗（v, e, i, 是否提前大修），按时间顺序。
// longest_path_w（T6）：提前大修弧只在启动窗口内存在，见文件末尾。
// 编译：g++ -O2 -shared -fPIC -o libprice.so price.cpp
#include <vector>
#include <cstring>
#include <limits>
extern "C" {

// w: ne × H（行主序）；nx: L × ne；nxR: ne。out_*: 长度至少 max_out。返回值岗个数（超过 max_out 返回 −1）。
int longest_path(int L, int H, int v0, int ne, const int* nx, const int* nxR, int reset,
                 const double* w, double* best_out,
                 int* out_v, int* out_e, int* out_i, int* out_r, int max_out) {
    const double NEG = -std::numeric_limits<double>::infinity();
    const int NN = L * H;
    std::vector<double> F(NN, NEG);
    // 前驱：节点编号、弧类型（0 等待，1 值岗，2 提前大修）、站 e
    std::vector<int> pn(NN, -1), pk(NN, -1), pe(NN, -1);
    double best = NEG; int bn = -1, bk = -1, be = -1;   // 进入汇点的最佳弧（尾节点、类型、站）
    F[0 * H + v0] = 0.0;
    for (int v = v0; v < H; v++) {
        for (int i = 0; i < L; i++) {
            const int node = i * H + v; const double f = F[node];
            if (f == NEG) continue;
            // 等待
            if (v + 1 < H) { int to = i * H + v + 1; if (f > F[to]) { F[to] = f; pn[to] = node; pk[to] = 0; pe[to] = -1; } }
            else if (f > best) { best = f; bn = node; bk = 0; be = -1; }
            for (int e = 0; e < ne; e++) {
                const double val = f + w[(size_t)e * H + v];
                // 正常值岗
                { int h = v + nx[i * ne + e], ni = (i + 1) % L;
                  if (h < H) { int to = ni * H + h; if (val > F[to]) { F[to] = val; pn[to] = node; pk[to] = 1; pe[to] = e; } }
                  else if (val > best) { best = val; bn = node; bk = 1; be = e; } }
                if (reset && i != L - 1) {
                    int h = v + nxR[e];
                    if (h < H) { int to = 0 * H + h; if (val > F[to]) { F[to] = val; pn[to] = node; pk[to] = 2; pe[to] = e; } }
                    else if (val > best) { best = val; bn = node; bk = 2; be = e; }
                }
            }
        }
    }
    *best_out = best;
    // 回溯：先收集到临时数组再倒序
    std::vector<int> tv, te, ti, tr;
    int node = bn, kind = bk, e = be;
    while (node >= 0) {
        if (kind >= 1) { tv.push_back(node % H); te.push_back(e); ti.push_back(node / H); tr.push_back(kind == 2); }
        int prev = pn[node]; int pkind = pk[node], pe2 = pe[node];
        if (prev < 0) break;
        kind = pkind; e = pe2; node = prev;
    }
    int n = (int)tv.size();
    if (n > max_out) return -1;
    for (int a = 0; a < n; a++) { out_v[a] = tv[n - 1 - a]; out_e[a] = te[n - 1 - a]; out_i[a] = ti[n - 1 - a]; out_r[a] = tr[n - 1 - a]; }
    return n;
}

// T6：与 longest_path 相同，只是提前大修弧只在值岗开始 D = v + ld[e] < W 时存在（启动窗口 [0, W)，对应 LP_fr^Rw）。
// ld: ne 个站的 lead。其余参数与输出同 longest_path。
int longest_path_w(int L, int H, int v0, int ne, const int* nx, const int* nxR, const int* ld, int W,
                   const double* w, double* best_out,
                   int* out_v, int* out_e, int* out_i, int* out_r, int max_out) {
    const double NEG = -std::numeric_limits<double>::infinity();
    const int NN = L * H;
    std::vector<double> F(NN, NEG);
    std::vector<int> pn(NN, -1), pk(NN, -1), pe(NN, -1);
    double best = NEG; int bn = -1, bk = -1, be = -1;
    F[0 * H + v0] = 0.0;
    for (int v = v0; v < H; v++) {
        for (int i = 0; i < L; i++) {
            const int node = i * H + v; const double f = F[node];
            if (f == NEG) continue;
            if (v + 1 < H) { int to = i * H + v + 1; if (f > F[to]) { F[to] = f; pn[to] = node; pk[to] = 0; pe[to] = -1; } }
            else if (f > best) { best = f; bn = node; bk = 0; be = -1; }
            for (int e = 0; e < ne; e++) {
                const double val = f + w[(size_t)e * H + v];
                { int h = v + nx[i * ne + e], ni = (i + 1) % L;
                  if (h < H) { int to = ni * H + h; if (val > F[to]) { F[to] = val; pn[to] = node; pk[to] = 1; pe[to] = e; } }
                  else if (val > best) { best = val; bn = node; bk = 1; be = e; } }
                if (i != L - 1 && v + ld[e] < W) {
                    int h = v + nxR[e];
                    if (h < H) { int to = 0 * H + h; if (val > F[to]) { F[to] = val; pn[to] = node; pk[to] = 2; pe[to] = e; } }
                    else if (val > best) { best = val; bn = node; bk = 2; be = e; }
                }
            }
        }
    }
    *best_out = best;
    std::vector<int> tv, te, ti, tr;
    int node = bn, kind = bk, e = be;
    while (node >= 0) {
        if (kind >= 1) { tv.push_back(node % H); te.push_back(e); ti.push_back(node / H); tr.push_back(kind == 2); }
        int prev = pn[node]; int pkind = pk[node], pe2 = pe[node];
        if (prev < 0) break;
        kind = pkind; e = pe2; node = prev;
    }
    int n = (int)tv.size();
    if (n > max_out) return -1;
    for (int a = 0; a < n; a++) { out_v[a] = tv[n - 1 - a]; out_e[a] = te[n - 1 - a]; out_i[a] = ti[n - 1 - a]; out_r[a] = tr[n - 1 - a]; }
    return n;
}
}
