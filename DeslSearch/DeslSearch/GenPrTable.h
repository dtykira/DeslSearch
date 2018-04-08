#pragma once
#include "types.h"
#include "GenPermTable.h"

#define MAX_OUTPUTDIFFS_NUMBER (8)

extern int PDT[SBOX_INPUTS_NUMBER][SBOX_OUTPUTS_NUMBER];

extern int PDT_Number;
extern int PDT_0_Number[SBOX_INPUTS_NUMBER][PR_NUMBER];
extern int PDT_0_Offset[SBOX_INPUTS_NUMBER][PR_NUMBER][2];
extern int PDT_1_Number[PR_NUMBER];
extern int PDT_1_Offset[PR_NUMBER][2];

extern int PDT_1_Non0Num[PR_NUMBER];
extern int PDT_1_Non0Val[PR_NUMBER][SBOX_INPUTS_NUMBER];

extern ALIGNED_TYPE_(u8,8) SPE[SBOX_NUMBER][SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER][SBOX_NUMBER];

extern u8 *sbox;
void Substitution(u8* output,u8 input);

void GenProb();
void GenPDT();
void Statistics();
void Storage();
void PDTStatistics();
void GenPrTable();