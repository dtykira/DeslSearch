#include "GenPrTable.h"

int PDT[SBOX_INPUTS_NUMBER][SBOX_OUTPUTS_NUMBER];
int PDT_Number;
int PDT_0_Number[SBOX_INPUTS_NUMBER][PR_NUMBER];
int PDT_0_Offset[SBOX_INPUTS_NUMBER][PR_NUMBER][2];
int PDT_1_Number[PR_NUMBER];
int PDT_1_Offset[PR_NUMBER][2];
int PDT_1_Non0Num[PR_NUMBER];
int PDT_1_Non0Val[PR_NUMBER][SBOX_INPUTS_NUMBER];

ALIGNED_TYPE_(u8,8) SPE[SBOX_NUMBER][SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER][SBOX_NUMBER];

#define PRINT_PDTTABLE 1
#define PRINT_PROB 1
#define PRINT_STATISTICS 1
#define PRINT_STORAGE 1
#define PRINT_PDTSTATISTICS 1

//替换函数S
u8 Sbox[4][16] ={
14,5,7,2,11,8,1,15,0,10,9,4,6,13,12,3,
5,0,8,15,14,3,2,12,11,7,6,9,13,4,1,10,
4,9,2,14,8,7,13,0,10,12,15,1,5,11,3,6,
9,6,15,5,3,8,4,11,7,1,12,2,0,14,10,13
};
u8 Substitution(u8 input){
	u8 x=input,x1,x2;
	x1=(x&0x1)|((x>>4)&0x2);
	x2=(x>>1)&0xf;
	return Sbox[x1][x2];
}

prType Prob[8];

//生成float型概率势数数组，概率从大到小
void GenProb(){
#if (PRINT_PROB)
	FILE *fp=fopen("Prob.txt","w");
#endif
	for(int i=0;i<=7;i++){
		Prob[i]=6.0-log((prType)((8-i)*2))/log(2.0);
#if (PRINT_PROB)
		fprintf(fp,"%f\t",Prob[i]);
#endif
	}
#if (PRINT_PROB)
	fclose(fp);
#endif
}

//生成概率分布表
void GenPDT(){
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		for(u16 j=0x0;j<SBOX_OUTPUTS_NUMBER;j++){
			PDT[i][j]=0x0;
		}
	}//全部初始化为0

	for(u16 id=0x0;id<SBOX_INPUTS_NUMBER;id++){
		for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
			u8 od=Substitution(i)^Substitution(i^id);
			PDT[id][od]++;
		}
	}

#if (PRINT_PDTTABLE)
	FILE *fp=fopen("PDT_Table.txt","w");
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"/*0x%02x*/{",i);
		for(u16 j=0x0;j<SBOX_OUTPUTS_NUMBER;j++){
			fprintf(fp,"%d,",PDT[i][j]);
		}
		fprintf(fp,"},\n");
	}
	fclose(fp);
#endif
}

//生成PDT统计信息，Storage()里用
void Statistics(){
	PDT_Number=0;
	memset(PDT_1_Number,0,PR_NUMBER*sizeof(int));
	memset(PDT_0_Number,0,PR_NUMBER*SBOX_INPUTS_NUMBER*sizeof(int));
	//初始化为0

	for(int i=0;i<SBOX_INPUTS_NUMBER;i++){
		for(int o=0;o<SBOX_OUTPUTS_NUMBER;o++){
			if(PDT[i][o]!=0&&PDT[i][o]!=64){
				PDT_Number++;
				PDT_1_Number[(8-PDT[i][o]/2)]++;
				PDT_0_Number[i][(8-PDT[i][o]/2)]++;
			}
		}
		PDT_0_Offset[i][0][0]=0;
		PDT_0_Offset[i][0][1]=PDT_0_Number[i][0];
		for(int p=1;p<PR_NUMBER;p++){
			PDT_0_Offset[i][p][0]=PDT_0_Offset[i][p-1][1];
			PDT_0_Offset[i][p][1]=PDT_0_Offset[i][p-1][1]+PDT_0_Number[i][p];
		}
	}
	PDT_1_Offset[0][0]=0;
	PDT_1_Offset[0][1]=PDT_1_Number[0];
	for(int p=1;p<PR_NUMBER;p++){
		PDT_1_Offset[p][0]=PDT_1_Offset[p-1][1];
		PDT_1_Offset[p][1]=PDT_1_Offset[p-1][1]+PDT_1_Number[p];
	}
}

void Storage(){
	int PDT_1_Count[PR_NUMBER];
	for(int p=0;p<PR_NUMBER;p++){
		PDT_1_Count[p]=PDT_1_Offset[p][0];
	}
	for(int i=0;i<SBOX_INPUTS_NUMBER;i++){
		int PDT_0_Count[PR_NUMBER];
		for(int p=0;p<PR_NUMBER;p++){
			PDT_0_Count[p]=PDT_0_Offset[i][p][0];
		}
		for(int o=0;o<SBOX_OUTPUTS_NUMBER;o++){
			if(PDT[i][o]!=0&&PDT[i][o]!=64){
				for(int si=0;si<SBOX_NUMBER;si++){
					memcpy(SPE[si][i][PDT_0_Count[8-PDT[i][o]/2]],PE_Table[si][o],SBOX_NUMBER);
				}
				PDT_1_Count[8-PDT[i][o]/2]++;
				PDT_0_Count[8-PDT[i][o]/2]++;
			}
		}
	}
}

void PDTStatistics(){
	for(int p=0;p<PR_NUMBER;p++){
		PDT_1_Non0Num[p]=0;
	}
	for(int p=0;p<PR_NUMBER;p++){
		for(u8 i=0;i<SBOX_INPUTS_NUMBER;i++){
			if(PDT_0_Number[i][p]!=0){
				PDT_1_Non0Val[p][PDT_1_Non0Num[p]]=i;
				PDT_1_Non0Num[0]++;
			}
		}
	}
}

void GenPrTable(){
	GenProb();
	GenPDT();
	Statistics();
	Storage();
	PDTStatistics();
}