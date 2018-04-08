#include "GenPermTable.h"

#define PRINT_PETABLE 1

ALIGNED_TYPE_(u8,8) PE_Table[SBOX_NUMBER][SBOX_INPUTS_NUMBER][SBOX_NUMBER];

void Gen_PE_Table(){
	ALIGNED_TYPE_(u8,8) Pid[SBOX_NUMBER],Pod[SBOX_NUMBER];
#if (PRINT_PETABLE)
	FILE *fp=fopen("PE_Table.txt","w");
#endif
	for(int s=0;s<SBOX_NUMBER;s++){
		memset(Pid,0,SBOX_NUMBER);
		memset(Pod,0,SBOX_NUMBER);
		for(u8 id=0;id<SBOX_INPUTS_NUMBER;id++){
			Pid[s]=id;
			PE(Pod,Pid);
			memcpy(PE_Table[s][id],Pod,sizeof(Pod));
#if (PRINT_PETABLE)
	fprintf(fp,"/*0x%02x*/{",id);
	for(u8 i=0;i<SBOX_NUMBER;i++){
		fprintf(fp,"0x%02x,",Pod[i]);
	}
	fprintf(fp,"},\n");
#endif
		}
#if (PRINT_PETABLE)
		fprintf(fp,"},\n�����ǵ�%d��S��\n",s);
#endif
	}
#if (PRINT_PETABLE)
	fprintf(fp,"}:\n\n");
	fclose(fp);
#endif
}