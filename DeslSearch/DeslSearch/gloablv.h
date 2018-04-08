#pragma once
#include "types.h"

extern bool firstRoundActive;

extern ALIGNED_TYPE_(u8,8) r_od[ROUND_N][SBOX_NUMBER];//�洢�ֲ��
extern si8 r_an[ROUND_N];//�洢r_od�Ļ�ԾS�и���
extern si8 r_ai[ROUND_N][SBOX_NUMBER];//�洢r_od�Ļ�ԾS�е��±�
extern prType r_pr[ROUND_N];

extern ALIGNED_TYPE_(u8,8) r_od_best[ROUND_N][SBOX_NUMBER];
extern si8 r_an_best[ROUND_N];
extern si8 r_ai_best[ROUND_N][SBOX_NUMBER];
extern prType r_pr_best[ROUND_N];

extern prType Bnc[ROUND_N];//Bnc[R]����R+1�ּ���ǰ��ѡweight
extern prType Bn[ROUND_N];//�������weight