#include "gloablv.h"

bool firstRoundActive;

ALIGNED_TYPE_(u8,8) r_od[ROUND_N][SBOX_NUMBER];
si8 r_an[ROUND_N];
si8 r_ai[ROUND_N][SBOX_NUMBER];
prType r_pr[ROUND_N];

ALIGNED_TYPE_(u8,8) r_od_best[ROUND_N][SBOX_NUMBER];
si8 r_an_best[ROUND_N];
si8 r_ai_best[ROUND_N][SBOX_NUMBER];
prType r_pr_best[ROUND_N];

prType Bnc[ROUND_N];
prType Bn[ROUND_N];