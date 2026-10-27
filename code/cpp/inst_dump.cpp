// 把 C++ 读入的实例逐项打印成 JSON，供 code/t1_check.py 与 Python 端逐项比对。
#include "mdrr_inst.h"
static string iv(const vector<int>& v) { string s = "["; for (size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + to_string(v[i]); return s + "]"; }
static string dv(double x) { char b[64]; snprintf(b, 64, "%.17g", x); return b; }
int main(int argc, char** argv) {
    for (int a = 1; a < argc; a++) {
        Inst I; I.load(argv[a]);
        printf("{\"file\":\"%s\",\"name\":\"%s\",\"nk\":%d,\"ns\":%d,\"nt\":%d,\"types\":[", argv[a], I.name.c_str(), I.nk, I.ns, I.nt);
        for (int k = 0; k < I.nk; k++) {
            const Type& t = I.T[k];
            printf("%s{\"d\":%d,\"L\":%d,\"v0\":%d,\"flex\":%d,\"elig\":%s,\"q\":%d,\"C\":%s,\"post\":%s,\"back\":%s,\"step\":[",
                   k ? "," : "", t.d, t.L, t.v0, t.flex, iv(t.elig).c_str(), t.q, iv(t.C).c_str(), dv(t.post_h).c_str(), dv(t.back_h).c_str());
            for (int i = 0; i < t.L; i++) { printf("%s[", i ? "," : ""); for (size_t e = 0; e < t.elig.size(); e++) printf("%s[%d,%d]", e ? "," : "", t.step[i][e].first, t.step[i][e].second); printf("]"); }
            printf("],\"pre\":%s}", iv(t.pre).c_str());
        }
        printf("],\"sites\":[");
        for (int s = 0; s < I.ns; s++) printf("%s[%s,%d,%d]", s ? "," : "", dv(I.S[s].delta).c_str(), I.S[s].lead, I.S[s].K);
        printf("],\"tools\":[");
        for (int j = 0; j < I.nt; j++) printf("%s[%s,%s,%s]", j ? "," : "", dv(I.TT[j].ins).c_str(), dv(I.TT[j].unin).c_str(), dv(I.TT[j].maint).c_str());
        printf("],\"N0\":%s,\"M0\":%s,\"nelig\":%s}\n", iv(I.N0).c_str(), iv(I.M0).c_str(), iv(I.nelig).c_str());
    }
}
