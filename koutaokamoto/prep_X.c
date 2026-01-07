#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRLEN 15
#define Q 3
#define SIGMA 10
#define QSIZE (SIGMA*SIGMA*SIGMA)

typedef struct {
    int *ids;
    int size;
    int cap;
} List;

// 位置付き3-gramインデックス
#define QS3 1000
#define POS3 (STRLEN - Q + 1)
// 位置付き4-gramインデックス
#define Q4 4
#define QS4 10000
List index3pos[POS3 * QS3];
#define POS4 (STRLEN - Q4 + 1)
// 提出用: アンカー位置はソース内固定（外部定義不要）
#undef SEG4_0
#undef SEG4_1
#undef SEG4_2
#undef SEG3_POS
#define SEG4_0 0
#define SEG4_1 4
#define SEG4_2 8
#define SEG3_POS (STRLEN - 3)

List index_tbl[QSIZE];
List index4[POS4 * QS4];
int db_size = 0;

static inline int encode3(const char *s) {
    int a = s[0]-'A'; if ((unsigned)a >= SIGMA) return -1;
    int b = s[1]-'A'; if ((unsigned)b >= SIGMA) return -1;
    int c = s[2]-'A'; if ((unsigned)c >= SIGMA) return -1;
    return a*100 + b*10 + c;
}

static inline int encode4(const char *s) {
    int a = s[0]-'A'; if ((unsigned)a >= SIGMA) return -1;
    int b = s[1]-'A'; if ((unsigned)b >= SIGMA) return -1;
    int c = s[2]-'A'; if ((unsigned)c >= SIGMA) return -1;
    int d = s[3]-'A'; if ((unsigned)d >= SIGMA) return -1;
    return a*1000 + b*100 + c*10 + d;
}


static inline void normalize15(char *s) {
    for (int i = 0; i < STRLEN; i++) {
        if (s[i] >= 'a' && s[i] <= 'j') s[i] = (char)(s[i] - 'a' + 'A');
    }
}

/* 2パス方式のため動的拡張は使用しない */

int main(int argc, char **argv) {
    if (argc != 2) return 1;

    FILE *fp = fopen(argv[1], "r");
    if (!fp) return 1;

    // 1st pass: カウントのみ
    int *cnt_tbl = calloc(QSIZE, sizeof(int));
    int *cnt3pos = calloc(POS3 * QS3, sizeof(int));
    int *cnt4    = calloc(POS4 * QS4, sizeof(int));
    if (!cnt_tbl || !cnt3pos || !cnt4) return 1;

    char buf[32];
    while (fgets(buf, sizeof(buf), fp)) {
        buf[strcspn(buf, "\n")] = 0;
        normalize15(buf);
        for (int i = 0; i <= STRLEN - Q; i++) {
            int k = encode3(&buf[i]);
            if (k >= 0) cnt_tbl[k]++;
        }
        // 使用する位置のみカウント（pos3はSEG3_POSのみ、pos4は3アンカーのみ）
        {
            int p = SEG3_POS; int k3p = encode3(&buf[p]); if (k3p >= 0) cnt3pos[p * QS3 + k3p]++;
        }
        {
            int p = SEG4_0; int k4 = encode4(&buf[p]); if (k4 >= 0) cnt4[p * QS4 + k4]++;
        }
        {
            int p = SEG4_1; int k4 = encode4(&buf[p]); if (k4 >= 0) cnt4[p * QS4 + k4]++;
        }
        {
            int p = SEG4_2; int k4 = encode4(&buf[p]); if (k4 >= 0) cnt4[p * QS4 + k4]++;
        }
        db_size++;
    }
    // 2nd passに備えて巻き戻し
    fseek(fp, 0, SEEK_SET);

    // 2nd pass用の連続メモリを確保してリストを割当
    long total_tbl = 0, total_3p = 0, total_4 = 0;
    for (int i = 0; i < QSIZE; i++) total_tbl += cnt_tbl[i];
    for (int i = 0; i < POS3 * QS3; i++) total_3p += cnt3pos[i];
    for (int i = 0; i < POS4 * QS4; i++) total_4  += cnt4[i];

    int *pool_tbl = total_tbl ? (int*)malloc(sizeof(int) * total_tbl) : NULL;
    int *pool_3p  = total_3p ? (int*)malloc(sizeof(int) * total_3p) : NULL;
    int *pool_4   = total_4  ? (int*)malloc(sizeof(int) * total_4 ) : NULL;
    if ((total_tbl && !pool_tbl) || (total_3p && !pool_3p) || (total_4 && !pool_4)) return 1;

    // プレフィックスサムでオフセット計算
    long off = 0; for (int i = 0; i < QSIZE; i++) { index_tbl[i].ids = pool_tbl + off; index_tbl[i].size = cnt_tbl[i]; index_tbl[i].cap = 0; off += cnt_tbl[i]; cnt_tbl[i] = 0; }
    off = 0; for (int i = 0; i < POS3*QS3; i++) { index3pos[i].ids = pool_3p + off; index3pos[i].size = cnt3pos[i]; index3pos[i].cap = 0; off += cnt3pos[i]; cnt3pos[i] = 0; }
    off = 0; for (int i = 0; i < POS4*QS4; i++) { index4[i].ids = pool_4 + off; index4[i].size = cnt4[i]; index4[i].cap = 0; off += cnt4[i]; cnt4[i] = 0; }

    // 出力: まずDBサイズ
    printf("%d\n", db_size);
    int id = 0;
    while (fgets(buf, sizeof(buf), fp)) {
        buf[strcspn(buf, "\n")] = 0;
        normalize15(buf);
        // DB行を即出力
        printf("%s\n", buf);
        // 各インデックスへ書込み
        for (int i = 0; i <= STRLEN - Q; i++) {
            int k = encode3(&buf[i]);
            if (k >= 0) index_tbl[k].ids[cnt_tbl[k]++] = id;
        }
        {
            int p = SEG3_POS; int k3p = encode3(&buf[p]); if (k3p >= 0) index3pos[p * QS3 + k3p].ids[cnt3pos[p * QS3 + k3p]++] = id;
        }
        {
            int p = SEG4_0; int k4 = encode4(&buf[p]); if (k4 >= 0) index4[p * QS4 + k4].ids[cnt4[p * QS4 + k4]++] = id;
        }
        {
            int p = SEG4_1; int k4 = encode4(&buf[p]); if (k4 >= 0) index4[p * QS4 + k4].ids[cnt4[p * QS4 + k4]++] = id;
        }
        {
            int p = SEG4_2; int k4 = encode4(&buf[p]); if (k4 >= 0) index4[p * QS4 + k4].ids[cnt4[p * QS4 + k4]++] = id;
        }
        id++;
    }
    fclose(fp);

    // グローバル3-gram索引
    for (int i = 0; i < QSIZE; i++) {
        printf("%d", index_tbl[i].size);
        for (int j = 0; j < index_tbl[i].size; j++) printf(" %d", index_tbl[i].ids[j]);
        printf("\n");
    }
    // 位置付き4-gram索引
    for (int p = 0; p < POS4; p++) {
        for (int k = 0; k < QS4; k++) {
            List *L = &index4[p * QS4 + k];
            printf("%d", L->size);
            for (int j = 0; j < L->size; j++) printf(" %d", L->ids[j]);
            printf("\n");
        }
    }
    // 位置付き3-gram索引
    for (int p = 0; p < POS3; p++) {
        for (int k = 0; k < QS3; k++) {
            List *L = &index3pos[p * QS3 + k];
            printf("%d", L->size);
            for (int j = 0; j < L->size; j++) printf(" %d", L->ids[j]);
            printf("\n");
        }
    }
    return 0;
}