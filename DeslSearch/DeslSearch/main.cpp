#include "types.h"
#include "GenPermTable.h"
#include "GenPrTable.h"
#include "search.h"
#include <stdint.h>
#define N0 (SBOX_NUMBER)
using namespace std;


int main(){
	Gen_PE_Table();
	GenPrTable();
	
	ALIGNED_TYPE_(u16,8) adult[8]={0x28, 0x02, 0, 0, 0, 0, 0, 0};
	__m128i *ta;
	ta=(__m128i *)adult;
	u16  *p = (u16 *) ta;
	for(int i=0;i<8;i++){
		printf("%x ",p[i]);
	}printf("\n");
	printf("bool:%d\n",isValid(*ta));

	__m128i zero=_mm_setzero_si128();
	__m128i *idp;
	__m128i tmp;
	idp=(__m128i *)SPE[0][1][0];

	__m128i tmp1;
	tmp1=_mm_set1_epi16(0x0030);
	printf("%d",isValid(*idp));
	
	/*FILE *fp=fopen("A0.txt","w");
	FILE *fq=fopen("C0.txt","w");
	si8 A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
	for(int M0=1;M0<=N0;M0++){
		for (si8 i=0; i<=(N0-M0); i++) A0[i] = 0;	for (si8 i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for (si8 i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1; T0[1] = 0; F0[N0] = N0 - M0 + 1; I0 = N0 - M0; L0 = N0;
		//初始化，前面的S盒不活跃，后面的M0个S盒活跃
		do
		{
			for(int i=1;i<=N0;i++){
				fprintf(fp,"%d ",A0[i]);
			}fprintf(fp,"\n");
			for(int i=1;i<=N0;i++){
				fprintf(fq,"%d ",C0[i]);
			}fprintf(fq,"\n");
				
			//生成下一个S盒活跃模式
			if (I0 == 0)
			{
				break;
			} 
			else
			{
				if (T0[I0] < 0) { if ((-T0[I0]) != (I0-1)){ T0[I0-1] = T0[I0]; } T0[I0] = I0-1; }
				if ( A0[I0]==0 )
				{
					X0 = I0; Y0 = F0[L0]; if (A0[I0-1] == 1){ F0[I0] = F0[I0 - 1]; } else { F0[I0] = I0; }
					if (F0[L0] == L0) { L0 = I0; I0 = T0[I0]; goto CHANGE0; }
					if (L0 == N0) { T0[F0[N0]] = -I0 - 1; T0[I0 + 1] = T0[I0]; I0 = F0[N0]; F0[N0] = F0[N0] + 1; goto CHANGE0; }
					T0[L0] = -I0-1; T0[I0+1] = T0[I0]; F0[L0] = F0[L0] + 1; I0 = L0; goto CHANGE0;
				}
				Y0 = I0;
				if (I0 != L0)
				{
					F0[L0] = X0 = F0[L0] - 1; F0[I0 - 1] = F0[I0];
					if (L0 == N0)
					{
						if (I0 == (F0[N0] - 1)) { I0 = T0[I0]; goto CHANGE0; }
						T0[F0[N0]-1] = -I0-1; T0[I0+1] = T0[I0]; I0 = F0[N0] - 1; goto CHANGE0;
					}
					T0[L0] = -I0 -1; T0[I0 + 1] = T0[I0]; I0 = L0; goto CHANGE0;
				}
				X0 = N0; F0[L0 - 1] = F0[L0]; F0[N0] = N0; L0 = N0;
				if (I0 == N0 - 1) { I0 = T0[N0 - 1]; goto CHANGE0; }
				T0[N0 - 1] = -I0 - 1; T0[I0 + 1] = T0[I0]; I0 = N0 - 1;
CHANGE0:
				A0[X0] = 1; A0[Y0] = 0; H0[X0] = Z0 = H0[Y0]; C0[Z0] = X0;
			}
		} while (true);
	}
	fclose(fp);
	fclose(fq);
	

	FILE *fw=fopen("a,m,j0.txt","w");
	// 这一组变量是为了以最小变动的次序遍历正向第一轮的输入差分的（在给定weight下至少有一个输入差分的输出差分）
	// 操作在其上的算法将会生成混合基系统的格雷Gray码（非二元的），新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	// 所需要关注的是 m0[], a0[], j0 
	// 其中，
	//     -- m0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的个数，初始化之后将不再变动
	//     -- a0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的下标，每次循环都会被改动，新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	//     -- j0 是整数，是每次循环中变动的那个元素的下标
	si8 m0[SBOX_NUMBER], a0[SBOX_NUMBER], d0[SBOX_NUMBER], f0[SBOX_NUMBER+1];
	si8 j0;
	si8 N1=3;

	j0 = 0;
	memset(a0, 0, sizeof(a0));
	for (si8 i=0; i<N1; i++) { d0[i] = 1; f0[i] = i;m0[i]=5-(i+1);} f0[N1] = N1;

	fprintf(fw,"m0:");
	for(int i=0;i<N1;i++){
		fprintf(fw,"%d ",m0[i]);
	}fprintf(fw,"\n");

	while (true)
	{	
		j0 = f0[0]; f0[0] = 0;

		fprintf(fw,"a0:");
		for(int i=0;i<N1;i++){
			fprintf(fw,"%d ",a0[i]);
		}fprintf(fw,"\tj0:%d\n",j0);

		if (j0==N1) break;
		else
		{
			a0[j0] = a0[j0] + d0[j0];
		}
		if ((a0[j0] == 0) || (a0[j0] == m0[j0]-1))
		{
			d0[j0] = -d0[j0];	f0[j0] = f0[j0+1];	f0[j0+1] = j0+1;
		}
	}
	fclose(fw);*/

	Round=2;
	Bnc[Round-1]=Bn[Round-2]+4;
	Round_1();

	system("pause");
	return 0;
}