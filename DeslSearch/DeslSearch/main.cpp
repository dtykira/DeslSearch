#include "types.h"
#include "GenPermTable.h"
#include "GenPrTable.h"

#define N0 (SBOX_NUMBER-1)

int main(){
	GenPrTable();

	// 这一组变量是为了以最小变动的次序遍历哪些S盒活跃的
	// 所需要关注的是 A0[] 和 C0[]，其它的是辅助的变量
	// 其中 A0[] 是 0 和 1 的序列，第 i 个是活跃的则 A[i] 为 1，否则为 0
	// 其中 C0[] 是 整数序列，第 i 个活跃S盒的下标可以通过查询 C0[i] 得到
	si8 A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
	for(int M0=1;M0<=SBOX_NUMBER;M0++){
		for (si8 i=0; i<=(N0-M0); i++) A0[i] = 0;	for (si8 i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for (si8 i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1; T0[1] = 0; F0[N0] = N0 - M0 + 1; I0 = N0 - M0; L0 = N0;

		do
		{
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
	system("pause");
	return 0;
}