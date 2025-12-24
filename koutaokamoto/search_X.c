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

typedef struct {
    int *ids;
    int size;
} List;

char **db;
List index_tbl[QSIZE];
int db_size;

int encode3(const char *s) {
    // 実データは 'A'..'J'（大文字）
    int a = s[0]-'A'; if ((unsigned)a >= SIGMA) return -1;
    int b = s[1]-'A'; if ((unsigned)b >= SIGMA) return -1;
    int c = s[2]-'A'; if ((unsigned)c >= SIGMA) return -1;
    return a*100 + b*10 + c;
}

/* クエリ用Peq（A..Jのみ）を構築 */
static inline void build_peq(const char *pattern, uint64_t peq[SIGMA]) {
    for (int c = 0; c < SIGMA; c++) peq[c] = 0;
    for (int i = 0; i < STRLEN; i++) {
        int idx = pattern[i] - 'A';
        if ((unsigned)idx < SIGMA) peq[idx] |= 1ULL << i;
    }
}

/* Myers法（編集距離 ≤3 判定）: peqはクエリから前計算済み */
static inline int myers_with_peq(const uint64_t *peq, const char *text) {
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
static inline int hamming_leq(const char *a, const char *b) {
    int diff = 0;
    for (int i = 0; i < STRLEN; i++) {
        diff += (a[i] != b[i]);
        if (diff > MAX_DIST) return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc != 3) return 1;

    FILE *fp = fopen(argv[2], "r");
    if (!fp) return 1;

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
    fclose(fp);

    FILE *fq = fopen(argv[1], "r");
    if (!fq) return 1;

    // スタンプ法で巨大memset回避
    uint32_t *seen = calloc(db_size, sizeof(uint32_t));
    uint8_t *cnt = (uint8_t*)malloc(sizeof(uint8_t) * db_size);
    int *cand = (int*)malloc(sizeof(int) * db_size);
    if (!seen || !cnt || !cand) return 1;
    uint32_t stamp = 1;

    while (fgets(buf, sizeof(buf), fq)) {
        buf[strcspn(buf, "\n")] = 0;
        int csz = 0;
        uint64_t peq[SIGMA];
        build_peq(buf, peq);
        // qグラムの重複をfreq集約し、ポスティングの小さい順に一度だけ走査
        const int POSN = STRLEN - Q + 1; // 13
        struct { int k; int sz; uint8_t freq; } pos[POSN];
        int u = 0;
        for (int i = 0; i < POSN; i++) {
            int k = encode3(&buf[i]);
            if (k < 0 || k >= QSIZE) continue;
            int found = -1;
            for (int t = 0; t < u; t++) if (pos[t].k == k) { found = t; break; }
            if (found >= 0) {
                if (pos[found].freq < 255) pos[found].freq++;
            } else {
                pos[u].k = k;
                pos[u].sz = index_tbl[k].size;
                pos[u].freq = 1;
                u++;
            }
        }
        // 挿入ソート（u<=13）
        for (int i = 1; i < u; i++) {
            int k = pos[i].k, sz = pos[i].sz; uint8_t fq = pos[i].freq;
            int j = i - 1;
            while (j >= 0 && pos[j].sz > sz) { pos[j+1] = pos[j]; j--; }
            pos[j+1].k = k; pos[j+1].sz = sz; pos[j+1].freq = fq;
        }
        // スタンプ法＋freq加算（しきい値到達で候補追加、重複追加は避ける）
        for (int pi = 0; pi < u; pi++) {
            int k = pos[pi].k; uint8_t fq = pos[pi].freq;
            int psz = index_tbl[k].size;
            int *ids = index_tbl[k].ids;
            for (int j = 0; j < psz; j++) {
                int id = ids[j];
                if (seen[id] != stamp) {
                    seen[id] = stamp;
                    uint8_t v = fq;
                    cnt[id] = v;
                    if (v >= MIN_MATCHES) cand[csz++] = id;
                } else {
                    uint8_t prev = cnt[id];
                    uint16_t nv = (uint16_t)prev + (uint16_t)fq;
                    if (nv > 255) nv = 255;
                    cnt[id] = (uint8_t)nv;
                    if (prev < MIN_MATCHES && nv >= MIN_MATCHES) cand[csz++] = id;
                }
            }
        }

        int found = 0;
        // 高カウント優先（降順）でMyersを適用し早期終了
        for (int need = POSN; need >= MIN_MATCHES && !found; need--) {
            for (int i = 0; i < csz; i++) {
                int id = cand[i];
                if (cnt[id] != need) continue;
                const char *s = db[id];
                if (hamming_leq(buf, s) || myers_with_peq(peq, s)) { found = 1; break; }
            }
        }
        putchar(found ? '1' : '0');
        // 次クエリ用にスタンプを進める（オーバーフロー時はゼロ化）
        if (++stamp == 0) {
            memset(seen, 0, sizeof(uint32_t) * db_size);
            stamp = 1;
        }
    }
    putchar('\n');
    return 0;
}