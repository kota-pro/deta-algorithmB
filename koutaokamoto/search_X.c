#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define STRLEN 15
#define Q 3
#define SIGMA 10
#define QSIZE (SIGMA*SIGMA*SIGMA)
#define MAX_DIST 3
// 編集距離kに対するqグラム共有下限: t = (m - q + 1) - k*q
#define MIN_MATCHES ((STRLEN - Q + 1) - (MAX_DIST * Q))
// 位置付きq-gram早期一致の位置バンド幅（±POS_BAND）
// 提出用: バンド幅はソース内で固定（外部定義不要）
#undef POS_BAND
#define POS_BAND 1
// 4-gramのアンカー位置をコンパイル時に指定（3点）
// 提出用: 4-gramアンカーは固定（0,4,8）
#undef SEG4_0
#undef SEG4_1
#undef SEG4_2
#define SEG4_0 0
#define SEG4_1 4
#define SEG4_2 8
// 3-gramのアンカー位置（既定は末尾ブロック12）
// 提出用: 3-gramアンカーは末尾（12）固定
#undef SEG3_POS
#define SEG3_POS (STRLEN - 3)
// 出力バッファサイズ（デフォルト64KB）
// 提出用: 出力バッファは32KB固定
#undef OUTBUF_SIZE
#define OUTBUF_SIZE 32768

typedef struct {
    int *ids;
    int size;
} List;

char **db;
List index_tbl[QSIZE];
// 位置付き4-gram索引
List index4[(STRLEN - 4 + 1) * 10000];
// 位置付き3-gram索引
List index3pos[(STRLEN - 3 + 1) * 1000];
// 位置付き5-gram索引（削除）
int db_size;

static inline int encode3(const char * restrict s) {
    // 実データは 'A'..'J'（大文字）
    int a = s[0]-'A'; if ((unsigned)a >= SIGMA) return -1;
    int b = s[1]-'A'; if ((unsigned)b >= SIGMA) return -1;
    int c = s[2]-'A'; if ((unsigned)c >= SIGMA) return -1;
    return a*100 + b*10 + c;
}

/* クエリ用Peq（A..Jのみ）を構築 */
static inline void build_peq(const char * restrict pattern, uint64_t peq[SIGMA]) {
    for (int c = 0; c < SIGMA; c++) peq[c] = 0;
    for (int i = 0; i < STRLEN; i++) {
        int idx = pattern[i] - 'A';
        if ((unsigned)idx < SIGMA) peq[idx] |= 1ULL << i;
    }
}

/* Myers法（編集距離 ≤3 判定）: peqはクエリから前計算済み */
static inline int myers_with_peq(const uint64_t * restrict peq, const char * restrict text) {
    uint64_t Pv = ~0ULL, Mv = 0;
    int score = STRLEN;
    for (int i = 0; i < STRLEN; i++) {
        int idx = text[i] - 'A';
        uint64_t Eq = (unsigned)idx < SIGMA ? peq[idx] : 0ULL;
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        if (Ph & (1ULL << (STRLEN - 1))) score++;
        if (Mh & (1ULL << (STRLEN - 1))) score--;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        if (score > MAX_DIST) return 0;
    }
    return score <= MAX_DIST;
}

/* ハミング距離が≤MAX_DISTなら即一致 */
static inline int hamming_leq(const char * restrict a, const char * restrict b) {
    int diff = 0;
    for (int i = 0; i < STRLEN; i++) {
        diff += (a[i] != b[i]);
        if (diff > MAX_DIST) return 0;
    }
    return 1;
}

/* 15文字の完全一致（memcpy-less, 未整列ロード想定） */
static inline int equal15(const char * restrict a, const char * restrict b) {
    const uint64_t *ua = (const uint64_t*)a; const uint64_t *ub = (const uint64_t*)b;
    if (*ua != *ub) return 0;
    const uint32_t *sa = (const uint32_t*)(a+8); const uint32_t *sb = (const uint32_t*)(b+8);
    if (*sa != *sb) return 0;
    const uint16_t *ta = (const uint16_t*)(a+12); const uint16_t *tb = (const uint16_t*)(b+12);
    if (*ta != *tb) return 0;
    return a[14] == b[14];
}

int main(int argc, char **argv) {
    if (argc != 3) return 1;

    FILE *fp = fopen(argv[2], "r");
    if (!fp) return 1;
    // 入出力のバッファ設定は標準のまま使用（環境依存の不具合回避）

    if (fscanf(fp, "%d", &db_size) != 1) return 1;
    int ch = fgetc(fp); (void)ch; // 改行を1つ消費
    db = malloc(sizeof(char*) * db_size);

    char buf[32];
    for (int i = 0; i < db_size; i++) {
        fgets(buf, sizeof(buf), fp);
        buf[strcspn(buf, "\n")] = 0;
        db[i] = strdup(buf);
    }

    for (int i = 0; i < QSIZE; i++) {
        fscanf(fp, "%d", &index_tbl[i].size);
        index_tbl[i].ids = malloc(sizeof(int) * index_tbl[i].size);
        for (int j = 0; j < index_tbl[i].size; j++) {
            fscanf(fp, "%d", &index_tbl[i].ids[j]);
        }
    }
    // 続けて位置付き4-gram索引（POS4=12, QS4=10000）を読み込む
    for (int p = 0; p < (STRLEN - 4 + 1); p++) {
        for (int k = 0; k < 10000; k++) {
            int idx = p * 10000 + k;
            fscanf(fp, "%d", &index4[idx].size);
            index4[idx].ids = malloc(sizeof(int) * index4[idx].size);
            for (int j = 0; j < index4[idx].size; j++) {
                fscanf(fp, "%d", &index4[idx].ids[j]);
            }
        }
    }
    // （5-gram読み込みは削除）
    // 続けて位置付き3-gram索引（POS3=13, QS3=1000）を読み込む
    for (int p = 0; p < (STRLEN - 3 + 1); p++) {
        for (int k = 0; k < 1000; k++) {
            int idx = p * 1000 + k;
            fscanf(fp, "%d", &index3pos[idx].size);
            index3pos[idx].ids = malloc(sizeof(int) * index3pos[idx].size);
            for (int j = 0; j < index3pos[idx].size; j++) {
                fscanf(fp, "%d", &index3pos[idx].ids[j]);
            }
        }
    }
    fclose(fp);

    FILE *fq = fopen(argv[1], "r");
    if (!fq) return 1;
    // 標準I/Oバッファはデフォルトのまま利用

    // スタンプ法で巨大memset回避（重複チェックのみ）
    uint16_t *seen_check = calloc(db_size, sizeof(uint16_t));
    if (!seen_check) return 1;
    uint16_t stamp = 1;

    // 安全な手動バッファ出力（適度なチャンク）
    char outbuf[OUTBUF_SIZE];
    size_t outlen = 0;
    while (fgets(buf, sizeof(buf), fq)) {
        buf[strcspn(buf, "\n")] = 0;
        // 入力は a..j も想定。A..J に正規化
        for (int i = 0; i < STRLEN; i++) {
            if (buf[i] >= 'a' && buf[i] <= 'j') buf[i] = (char)(buf[i] - 'a' + 'A');
        }

        uint64_t peq[SIGMA];
        build_peq(buf, peq);

        // 位置付き4-gram/3-gramの早期一致チェック（±2の位置バンド）
        int found = 0;
        // 4-gram（pos: SEG4_0, SEG4_1, SEG4_2）
        const int POS4 = STRLEN - 4 + 1; // 12
        const int seg_pos4[3] = {SEG4_0, SEG4_1, SEG4_2};
        for (int s = 0; s < 3 && !found; s++) {
            int pos = seg_pos4[s];
            int a = buf[pos]-'A'; int b = buf[pos+1]-'A'; int c = buf[pos+2]-'A'; int d = buf[pos+3]-'A';
            if ((unsigned)a>=SIGMA||(unsigned)b>=SIGMA||(unsigned)c>=SIGMA||(unsigned)d>=SIGMA) continue;
            int key4 = a*1000 + b*100 + c*10 + d;
            int lo = pos - POS_BAND; if (lo < 0) lo = 0;
            int hi = pos + POS_BAND; if (hi >= POS4) hi = POS4 - 1;
            for (int p = lo; p <= hi && !found; p++) {
                int idx = p * 10000 + key4;
                int psz = index4[idx].size;
                int *ids = index4[idx].ids;
                for (int j = 0; j < psz; j++) {
                    int id = ids[j];
                    if (__builtin_expect(seen_check[id] != stamp, 1)) {
                        seen_check[id] = stamp;
                        const char *sdb = db[id];
                        if (equal15(buf, sdb) || hamming_leq(buf, sdb) || myers_with_peq(peq, sdb)) { found = 1; break; }
                    }
                }
            }
        }
        // 3-gram（アンカー位置 SEG3_POS）
        if (!found) {
            int pos3 = SEG3_POS; // 12既定
            int a = buf[pos3]-'A'; int b = buf[pos3+1]-'A'; int c = buf[pos3+2]-'A';
            if ((unsigned)a<SIGMA&&(unsigned)b<SIGMA&&(unsigned)c<SIGMA) {
                int key3 = a*100 + b*10 + c;
                int lo = pos3 - POS_BAND; if (lo < 0) lo = 0;
                int hi = pos3 + POS_BAND; int POS3N = STRLEN - 3 + 1; if (hi >= POS3N) hi = POS3N - 1;
                for (int p = lo; p <= hi && !found; p++) {
                    int idx = p * 1000 + key3;
                    int psz = index3pos[idx].size;
                    int *ids = index3pos[idx].ids;
                    for (int j = 0; j < psz; j++) {
                        int id = ids[j];
                        if (__builtin_expect(seen_check[id] != stamp, 1)) {
                            seen_check[id] = stamp;
                            const char *sdb = db[id];
                            if (equal15(buf, sdb) || hamming_leq(buf, sdb) || myers_with_peq(peq, sdb)) { found = 1; break; }
                        }
                    }
                }
            }
        }
        if (found) {
            outbuf[outlen++] = '1'; if (outlen == sizeof(outbuf)) { fwrite(outbuf, 1, outlen, stdout); outlen = 0; }
            if (++stamp == 0) { memset(seen_check, 0, sizeof(uint16_t) * db_size); stamp = 1; }
            continue;
        }

        // q=3フィルタ：最小2バケットの交差＋3本目で絞り、残りは二分探索
        const int POSN = STRLEN - Q + 1; // 13
        struct { int k; int sz; const int *ids; } pos[POSN];
        int u = 0;
        for (int i = 0; i < POSN; i++) {
            int k = encode3(&buf[i]);
            if (k < 0 || k >= QSIZE) continue;
            pos[u].k = k;
            pos[u].sz = index_tbl[k].size;
            pos[u].ids = index_tbl[k].ids;
            u++;
        }
        // 3最小のみ選択（全体ソートは避ける）
        int m1=-1,m2=-1,m3=-1;
        for (int i = 0; i < u; i++) {
            if (m1<0 || pos[i].sz < pos[m1].sz) { m3 = m2; m2 = m1; m1 = i; }
            else if (m2<0 || pos[i].sz < pos[m2].sz) { m3 = m2; m2 = i; }
            else if (m3<0 || pos[i].sz < pos[m3].sz) { m3 = i; }
        }

            int found2 = 0;
            if (u >= 3) {
                // 2リスト交差→3リスト目で絞る（m1,m2,m3）
                const int *A = pos[m1].ids; int na = pos[m1].sz;
                const int *B = pos[m2].ids; int nb = pos[m2].sz;
                int ia = 0, ib = 0, inter_sz = 0;
                static int *inter = NULL; static int inter_cap = 0;
                int need = (na < nb ? na : nb);
                if (need > inter_cap) {
                    int *tmp = (int*)realloc(inter, sizeof(int) * need);
                    if (!tmp) return 1;
                    inter = tmp; inter_cap = need;
                }
                while (ia < na && ib < nb) {
                    int va = A[ia], vb = B[ib];
                    if (va == vb) { inter[inter_sz++] = va; ia++; ib++; }
                    else if (va < vb) ia++; else ib++;
                }
                // 3本目と交差（2ポインタかバイナリ検索）
                const int *C = pos[m3].ids; int nc = pos[m3].sz;
                int ic = 0; // Cのポインタ

                for (int t = 0; t < inter_sz && !found2; t++) {
                    int v = inter[t];
                    while (ic < nc && C[ic] < v) ic++;
                    if (ic < nc && C[ic] == v) {
                        // すでに3ヒット。残りでしきい値到達を確認
                        int cntv = 3;
                        for (int pi = 0; pi < u && cntv < MIN_MATCHES; pi++) {
                            if (pi==m1||pi==m2||pi==m3) continue;
                            int l = 0, r = pos[pi].sz - 1; const int *arr = pos[pi].ids; int ok = 0;
                            while (l <= r) {
                                int m = (l + r) >> 1; int x = arr[m];
                                if (x == v) { ok = 1; break; }
                                if (x < v) l = m + 1; else r = m - 1;
                            }
                            if (ok) cntv++;
                        }
                        if (cntv >= MIN_MATCHES) {
                            const char *sdb = db[v];
                            if (equal15(buf, sdb) || hamming_leq(buf, sdb) || myers_with_peq(peq, sdb)) { found2 = 1; break; }
                        }
                    }
                }
            } else if (u >= 2) {
                // 2リスト交差のみ
                // 最小2本を使用
                int m1=-1,m2=-1; for (int i=0;i<u;i++){ if(m1<0||pos[i].sz<pos[m1].sz){ m2=m1; m1=i; } else if(m2<0||pos[i].sz<pos[m2].sz){ m2=i; }}
                const int *A = pos[m1].ids; int na = pos[m1].sz;
                const int *B = pos[m2].ids; int nb = pos[m2].sz;
                int ia = 0, ib = 0;
                while (ia < na && ib < nb && !found2) {
                    int va = A[ia], vb = B[ib];
                    if (va == vb) {
                        int cntv = 2;
                        for (int pi = 0; pi < u && cntv < MIN_MATCHES; pi++) {
                            if (pi==m1||pi==m2) continue;
                            int l = 0, r = pos[pi].sz - 1; const int *arr = pos[pi].ids; int ok = 0;
                            while (l <= r) {
                                int m = (l + r) >> 1; int x = arr[m];
                                if (x == va) { ok = 1; break; }
                                if (x < va) l = m + 1; else r = m - 1;
                            }
                            if (ok) cntv++;
                        }
                        if (cntv >= MIN_MATCHES) {
                            const char *sdb = db[va];
                            if (equal15(buf, sdb) || hamming_leq(buf, sdb) || myers_with_peq(peq, sdb)) { found2 = 1; break; }
                        }
                        ia++; ib++;
                    } else if (va < vb) {
                        ia++;
                    } else {
                        ib++;
                    }
                }
            }
            outbuf[outlen++] = (char)(found2 ? '1' : '0');
            if (outlen == sizeof(outbuf)) { fwrite(outbuf, 1, outlen, stdout); outlen = 0; }
        // 次クエリ用にスタンプを進める（オーバーフロー時はゼロ化）
        if (++stamp == 0) {
            memset(seen_check, 0, sizeof(uint16_t) * db_size);
              stamp = 1;
        }
    }
    if (outlen) fwrite(outbuf, 1, outlen, stdout);
    fputc('\n', stdout);
    return 0;
}