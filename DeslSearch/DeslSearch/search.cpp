#include "search.h"

#define N0 (SBOX_NUMBER-1)
#define Ci(i) C0[i]//(C0[i+1]-1)
FILE *fp_trails;

void print128(__m128i tmp){
	u16  *p = (u16 *) &tmp;
	for(int i=0;i<8;i++){
		printf("%02x ",p[i]);
	}printf("\n");
}

void fprint128(__m128i tmp){
	u16  *p = (u16 *) &tmp;
	for(int i=0;i<8;i++){
		fprintf(fp_trails,"%02x ",p[i]);
	}fprintf(fp_trails,"\n");
}

__m128i rshift(__m128i x,int i){
	__m128i tmp1,tmp2;
	switch(i){
	case 1:
		tmp1=_mm_slli_si128(x,2);//右移i个s盒
		tmp2=_mm_srli_si128(x,14);//左移8-i个S盒
		break;
	case 2:
		tmp1=_mm_slli_si128(x,4);//右移i个s盒
		tmp2=_mm_srli_si128(x,12);//左移8-i个S盒
		break;
	case 3:
		tmp1=_mm_slli_si128(x,6);//右移i个s盒
		tmp2=_mm_srli_si128(x,10);//左移8-i个S盒
		break;
	case 4:
		tmp1=_mm_slli_si128(x,8);//右移i个s盒
		tmp2=_mm_srli_si128(x,8);//左移8-i个S盒
		break;
	case 5:
		tmp1=_mm_slli_si128(x,10);//右移i个s盒
		tmp2=_mm_srli_si128(x,6);//左移8-i个S盒
		break;
	case 6:
		tmp1=_mm_slli_si128(x,12);//右移i个s盒
		tmp2=_mm_srli_si128(x,4);//左移8-i个S盒
		break;
	case 7:
		tmp1=_mm_slli_si128(x,14);//右移i个s盒
		tmp2=_mm_srli_si128(x,2);//左移8-i个S盒
		break;
	}
	return _mm_xor_si128(tmp1,tmp2);
}

bool isDup(__m128i tmp){
	u16 *p = (u16 *) &tmp;
	return p[0]!=0;
}

bool isValid(__m128i tmp){
	__m128i mask;
	mask=_mm_set1_epi16(0x0003);
	
	__m128i tmp0,tmp1,tmp2,tmp3,tmpCmp,tmpRes,Result;

	//返回一个__m128i的寄存器，将寄存器_A中的8个16bit整数按照_Count进行相同的逻辑右移，
	//移位填充值为0,r0=srl(_A0, _Count), r1=srl(_A1, _Count), ... r7=srl(_A7, _Count), 
	//shifting in zeros
	tmp0=_mm_srli_epi16(tmp, 4);// ||**** → 0000||

	//返回一个__m128i的寄存器，r=srl(_A, _Imm * 8),   _Imm must be an immediate,  
	//shifting in zeros
	//print128(tmp);
	//print128(tmp0);
	//对于这里si128的字节移位，是方向反过来的。
	tmp1=_mm_srli_si128(tmp0, 2);//print128(tmp1);
	tmp2=_mm_slli_si128(tmp0, 14);//print128(tmp2);
	tmp3=_mm_xor_si128(tmp1,tmp2);//print128(tmp3);
	
	tmpCmp=_mm_xor_si128(tmp,tmp3);//print128(tmpCmp);
	tmpRes=_mm_and_si128(tmpCmp,mask);
	
	u64 *p = (u64 *) &tmpRes;
	return (p[0]==0&&p[1]==0);
}

__m128i transform(u16 idv,si8 si){
	ALIGNED_TYPE_(u16,8) t[8];
	memset(t,0,16);
	t[si]=idv;
	__m128i *tmp;
	tmp=(__m128i *)t;
	return *tmp;
}



void setAndPrint(){
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[0][si]);
	}fprintf(fp_trails,"\t%f\n",r_pr[0]);

	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[1][si]);
	}fprintf(fp_trails,"\t");
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_r[1][si]);
	}fprintf(fp_trails,"\t%f\n",r_pr[1]);

	for(int r=2;r<Round-1;r++){
		for(int si=0;si<SBOX_NUMBER;si++){
			fprintf(fp_trails,"%02x ",r_od_l[r][si]);
		}fprintf(fp_trails,"\t");
		for(int si=0;si<SBOX_NUMBER;si++){
			fprintf(fp_trails,"%02x ",r_od_r[r][si]);
		}
		fprintf(fp_trails,"\t%f",r_pr[r]);
		fprintf(fp_trails,"\t%d: ",r_an[r]);
		for(int si=0;si<r_an[r];si++){
			fprintf(fp_trails,"%d ",r_ai[r][si]);
		}fprintf(fp_trails,"\n");
	}
	
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[Round-1][si]);
	}
	fprintf(fp_trails,"\t%f",r_pr[Round-1]);
	fprintf(fp_trails,"\t%d: ",r_an[Round-1]);
	for(int si=0;si<r_an[Round-1];si++){
		fprintf(fp_trails,"%d ",r_ai[Round-1][si]);
	}fprintf(fp_trails,"\n");
	
	fprintf(fp_trails,"---------------\n");
}

void getInfo(int r,__m128i tmp){
	__m128i a;
	__m128i b;
	__m128i c;
	__m128i d;
	int lm;
	int la;

	//获得r_an[r]
	a=_mm_setzero_si128();
	a=_mm_cmpgt_epi16(tmp,a);//tmp非零的16比特置为0xffff
	lm=_mm_movemask_epi8(a);
	la=_mm_popcnt_u32(lm)/2;
	r_an[r]=la;

	//获得r_ai[r][]
	c=_mm_load_si128((__m128i *)(&(W16v[lm][0])));
	_mm_storeu_si128((__m128i *)(r_ai[r]),c);
}

void Round_last(){
	__m128i *dp1;
	__m128i *dp2;
	__m128i *idp;
	__m128i tmp0;
	prType pr=0;

	dp1=(__m128i *)(r_od_l[Round-3]);
	dp2=(__m128i *)(r_od_r[Round-2]);
	idp=(__m128i *)(r_od_l[Round-1]);
	tmp0=_mm_xor_si128(*dp1,*dp2);
	_mm_store_si128(idp,tmp0);
	getInfo(Round-1,*idp);
	for(int si=0;si<r_an[Round-1];si++){
		pr+=PDT_MaxProb[r_od_l[Round-1][r_ai[Round-1][si]]];
	}
	if(pr+r_pr[Round-2]<(Bnc[Round-1]+1e-10)){
		r_pr[Round-1]=pr+r_pr[Round-2];
		setAndPrint();
	}
}

void Round_2(){
	clock_t start,end;
	start=clock();
	__m128i *idp;
	__m128i *odp;
	idp=(__m128i *)(r_od_l[1]);//idp指向r_od_l
	odp=(__m128i *)(r_od_r[1]);//odp指向r_od_r
	__m128i tmp0;//存放输入差分的局部变量
	__m128i tmp1;//存放输出差分的局部变量
	__m128i tmp2[8];
	prType pr[SBOX_NUMBER+1]={0};//计算概率用
	

	//0个S盒
	tmp0=_mm_setzero_si128();
	_mm_store_si128(idp,tmp0);
	_mm_store_si128(odp,tmp0);
	r_pr[1]=r_pr[0];
	if(firstRoundActive){Round_last();}
	//1个S盒
	for(int i=0;i<SBOX_INPUTS_NUMBER-1;i++){
		pr[1]=PDT_MaxProb[WtiForTravel[i]];
		if((pr[1]+Bn[Round-3]+r_pr[0])>=(Bnc[Round-1]+1e-10)){break;}
		tmp0=transform(WtiForTravel[i],0);
		if(!isValid(tmp0)){continue;}
		//_mm_store_si128(idp,tmp0);

		for(int j=0;j<WtoForTravelNumber[WtiForTravel[i]];j++){
			_mm_store_si128(idp,tmp0);
			tmp1=_mm_load_si128((__m128i *)(SPE[0][WtiForTravel[i]][j]));
			pr[1]=WtoForTravelProb[WtiForTravel[i]][j];
			if((pr[1]+Bn[Round-3]+r_pr[0])>=(Bnc[Round-1]+1e-10)){break;}

			_mm_store_si128(odp,tmp1);
			r_pr[1]=pr[1]+r_pr[0];
			Round_last();
			for(int k=1;k<=7;k++){
				_mm_store_si128(idp,rshift(tmp0,k));
				tmp1=_mm_load_si128((__m128i *)(SPE[k][WtiForTravel[i]][j]));
				_mm_store_si128(odp,tmp1);
				Round_last();
			}
		}
	}
	//活跃S盒两个及以上
	
	// 这一组变量是为了以最小变动的次序遍历哪些S盒活跃的
	// 所需要关注的是 A0[] 和 C0[]，其它的是辅助的变量
	// 其中 A0[] 是 0 和 1 的序列，第 i 个是活跃的则 A[i] 为 1，否则为 0
	// 其中 C0[] 是 整数序列，第 i 个活跃S盒的下标可以通过查询 C0[i] 得到
	si8 A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
	C0[0]=0;
	A0[0]=1;
	// 这一组变量是为了以最小变动的次序遍历正向第一轮的输入差分的（在给定weight下至少有一个输入差分的输出差分）
	// 操作在其上的算法将会生成混合基系统的格雷Gray码（非二元的），新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	// 所需要关注的是 m0[], a0[], j0 
	// 其中，
	//     -- m0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的个数，初始化之后将不再变动
	//     -- a0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的下标，每次循环都会被改动，新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	//     -- j0 是整数，是每次循环中变动的那个元素的下标
	si8 m0[SBOX_NUMBER], a0[SBOX_NUMBER], d0[SBOX_NUMBER], f0[SBOX_NUMBER+1];
	si8 j0=0;
	
	si8 m1[SBOX_NUMBER], a1[SBOX_NUMBER], d1[SBOX_NUMBER], f1[SBOX_NUMBER+1];
	si8 j1=0;

	si8 N1;
	//存放输入差分的局部变量
	si8 idv[SBOX_NUMBER];

	for(int M0=1;M0<=7;M0++){//活跃S盒两个及以上.M0<=N0
		if(((M0+1)*WMIN_S+Bn[Round-2])>=(Bnc[Round-1]+1e-10)){break;}
		N1=M0+1;
		for (si8 i=0; i<=(N0-M0); i++) A0[i] = 0;	for (si8 i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for (si8 i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1; T0[1] = 0; F0[N0] = N0 - M0 + 1; I0 = N0 - M0; L0 = N0;
		//初始化，前面的S盒不活跃，后面的M0个S盒活跃
		do
		{
			//生成第一个输入差分
			tmp0=_mm_setzero_si128();
			tmp1=_mm_setzero_si128();
			for(int l=1;l<=7;l++){
				tmp2[l]=_mm_setzero_si128();
			}

			for(si8 i=0;i<N1;i++){
				m0[i]=SBOX_INPUTS_NUMBER-1;
				idv[i]=WtiForTravel[0];
				tmp0=_mm_xor_si128(tmp0,transform(idv[i],Ci(i)));

				//对第一个输入差分，生成第一个输出差分
				tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(i)][idv[i]][0]));
				for(int l=1;l<=7;l++){
					tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(i)+l)%8][idv[i]][0]));
				}
				m1[i]=WtoForTravelNumber[idv[i]];

				pr[i+1]=pr[i]+PDT_MaxProb[idv[i]];
			}
			if(isValid(tmp0)&&(pr[N1]+Bn[Round-3]+r_pr[0])<(Bnc[Round-1]+1e-10)){
				r_pr[1]=pr[N1]+r_pr[0];
				_mm_store_si128(idp,tmp0);

				//对第一个输入差分，遍历输出差分
				j1 = 0;
				memset(a1, 0, sizeof(a1));
				for (si8 i=0; i<N1; i++) { d1[i] = 1; f1[i] = i;} f1[N1] = N1;
				while (true)
				{	
					j1 = f1[0]; f1[0] = 0;
					if (j1==N1) break;

					tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(j1)][idv[j1]][a1[j1]]));
					for(int l=1;l<=7;l++){
						tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(j1)+l)%8][idv[j1]][a1[j1]]));
					}
					pr[N1]-=WtoForTravelProb[idv[j1]][a1[j1]];

					a1[j1] = a1[j1] + d1[j1];

					if ((a1[j1] == 0) || (a1[j1] == m1[j1]-1))
					{
						d1[j1] = -d1[j1];	f1[j1] = f1[j1+1];	f1[j1+1] = j1+1;
					}

					tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(j1)][idv[j1]][a1[j1]]));
					for(int l=1;l<=7;l++){
						tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(j1)+l)%8][idv[j1]][a1[j1]]));
					}
					pr[N1]+=WtoForTravelProb[idv[j1]][a1[j1]];

					if(((pr[N1]+Bn[Round-3]+r_pr[0])<(Bnc[Round-1]+1e-10))){
						r_pr[1]=pr[N1]+r_pr[0];
						_mm_store_si128(idp,tmp0);
						_mm_store_si128(odp,tmp1);
						Round_last();
						for(int l=1;l<=7;l++){
							if(isDup(rshift(tmp0,l))){continue;}
							_mm_store_si128(idp,rshift(tmp0,l));
							_mm_store_si128(odp,tmp2[l]);
							Round_last();
						}
					}
				}
			}

			//初始化数组
			j0 = 0;
			memset(a0, 0, sizeof(a0));
			for (si8 i=0; i<N1; i++) { d0[i] = 1; f0[i] = i;} f0[N1] = N1;
			while (true)
			{	
				j0 = f0[0]; f0[0] = 0;
				if (j0==N1) break;
				else
				{
					idv[j0]=WtiForTravel[a0[j0]];
					tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));
					pr[N1]-=PDT_MaxProb[idv[j0]];
					a0[j0] = a0[j0] + d0[j0];
				}
				if ((a0[j0] == 0) || (a0[j0] == m0[j0]-1))
				{
					d0[j0] = -d0[j0];	f0[j0] = f0[j0+1];	f0[j0+1] = j0+1;
				}
				idv[j0]=WtiForTravel[a0[j0]];
				tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));

				//对每一个有效的输入差分，生成第一个输出差分
				tmp1=_mm_setzero_si128();
				for(int l=1;l<=7;l++){
					tmp2[l]=_mm_setzero_si128();
				}
				for(si8 i=0;i<N1;i++){
					tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(i)][idv[i]][0]));
					for(int l=1;l<=7;l++){
						tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(i)+l)%8][idv[i]][0]));
					}
					m1[i]=WtoForTravelNumber[idv[i]];
				}
				

				pr[N1]+=PDT_MaxProb[idv[j0]];

				if(isValid(tmp0)&&((pr[N1]+Bn[Round-3]+r_pr[0])<(Bnc[Round-1]+1e-10))){

					//对每一个有效的输入差分，遍历输出差分
					j1 = 0;
					memset(a1, 0, sizeof(a1));
					for (si8 i=0; i<N1; i++) { d1[i] = 1; f1[i] = i;} f1[N1] = N1;
					while (true)
					{	
						j1 = f1[0]; f1[0] = 0;
						if (j1==N1) break;

						tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(j1)][idv[j1]][a1[j1]]));
						for(int l=1;l<=7;l++){
							tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(j1)+l)%8][idv[j1]][a1[j1]]));
						}
						pr[N1]-=WtoForTravelProb[idv[j1]][a1[j1]];

						a1[j1] = a1[j1] + d1[j1];

						if ((a1[j1] == 0) || (a1[j1] == m1[j1]-1))
						{
							d1[j1] = -d1[j1];	f1[j1] = f1[j1+1];	f1[j1+1] = j1+1;
						}

						tmp1=_mm_xor_si128(tmp1,*(__m128i *)(SPE[Ci(j1)][idv[j1]][a1[j1]]));
						for(int l=1;l<=7;l++){
							tmp2[l]=_mm_xor_si128(tmp2[l],*(__m128i *)(SPE[(Ci(j1)+l)%8][idv[j1]][a1[j1]]));
						}
						pr[N1]+=WtoForTravelProb[idv[j1]][a1[j1]];

						if(((pr[N1]+Bn[Round-3]+r_pr[0])<(Bnc[Round-1]+1e-10))){
							r_pr[1]=pr[N1]+r_pr[0];
							_mm_store_si128(idp,tmp0);
							_mm_store_si128(odp,tmp1);
							Round_last();
							for(int l=1;l<=7;l++){
								if(isDup(rshift(tmp0,l))){continue;}
								_mm_store_si128(idp,rshift(tmp0,l));
								_mm_store_si128(odp,tmp2[l]);
								Round_last();
							}
						}
					}
				}
			}



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
	end = clock();
	timeForRound2+=((double)(end-start)/CLK_TCK);
}

void Round_1(){
	fp_trails=fopen("trails.txt","w");
	__m128i *idp;
	idp=(__m128i *)(r_od_l[0]);//idp指向r_od_l
	__m128i tmp0;//存放输入差分的局部变量
	__m128i tmp;
	prType pr[SBOX_NUMBER+1]={0};//计算概率用
	

	//0个S盒
	//memset(r_od_l[0],0,SBOX_NUMBER*sizeof(u16));
	*idp=_mm_setzero_si128();
	firstRoundActive=false;
	r_pr[0]=0;
	Round_2();

	//1个S盒
	firstRoundActive=true;
	for(int i=0;i<SBOX_INPUTS_NUMBER-1;i++){
		pr[1]=PDT_MaxProb[WtiForTravel[i]];
		if((pr[1]+Bn[Round-2])>=(Bnc[Round-1]+1e-10)){break;}
		tmp0=transform(WtiForTravel[i],0);
		if(!isValid(tmp0)){continue;}
		_mm_store_si128(idp,tmp0);
		r_pr[0]=pr[1];
		Round_2();
		for(int i=1;i<=7;i++){
			_mm_store_si128(idp,rshift(tmp0,i));
			Round_2();
		}
	}
	//活跃S盒两个及以上
	
	// 这一组变量是为了以最小变动的次序遍历哪些S盒活跃的
	// 所需要关注的是 A0[] 和 C0[]，其它的是辅助的变量
	// 其中 A0[] 是 0 和 1 的序列，第 i 个是活跃的则 A[i] 为 1，否则为 0
	// 其中 C0[] 是 整数序列，第 i 个活跃S盒的下标可以通过查询 C0[i] 得到
	si8 A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
	C0[0]=0;
	A0[0]=1;
	// 这一组变量是为了以最小变动的次序遍历正向第一轮的输入差分的（在给定weight下至少有一个输入差分的输出差分）
	// 操作在其上的算法将会生成混合基系统的格雷Gray码（非二元的），新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	// 所需要关注的是 m0[], a0[], j0 
	// 其中，
	//     -- m0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的个数，初始化之后将不再变动
	//     -- a0[] 是整数序列，是每个活跃S盒在给定weight下至少有一个输入差分的输出差分的下标，每次循环都会被改动，新的一组 a0[] 在前一组 a0[] 的基础上只变动一个元素
	//     -- j0 是整数，是每次循环中变动的那个元素的下标
	si8 m0[SBOX_NUMBER], a0[SBOX_NUMBER], d0[SBOX_NUMBER], f0[SBOX_NUMBER+1];
	si8 j0;
	j0 = 0;
	si8 N1;
	//存放输入差分的局部变量
	si8 idv[SBOX_NUMBER];

	for(int M0=1;M0<=7;M0++){//活跃S盒两个及以上.M0<=N0
		if(((M0+1)*WMIN_S+Bn[Round-2])>=(Bnc[Round-1]+1e-10)){break;}
		N1=M0+1;
		for (si8 i=0; i<=(N0-M0); i++) A0[i] = 0;	for (si8 i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for (si8 i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1; T0[1] = 0; F0[N0] = N0 - M0 + 1; I0 = N0 - M0; L0 = N0;
		//初始化，前面的S盒不活跃，后面的M0个S盒活跃
		do
		{
			//生成第一个输入差分
			tmp0=_mm_setzero_si128();
			for(si8 i=0;i<N1;i++){
				m0[i]=SBOX_INPUTS_NUMBER-1;
				idv[i]=WtiForTravel[0];
				tmp0=_mm_xor_si128(tmp0,transform(idv[i],Ci(i)));
				pr[i+1]=pr[i]+PDT_MaxProb[idv[i]];
			}
			if(isValid(tmp0)&&(pr[N1]+Bn[Round-2])<(Bnc[Round-1]+1e-10)){
				r_pr[0]=pr[N1];
				_mm_store_si128(idp,tmp0);
				Round_2();
				for(int i=1;i<=7;i++){
					if(isDup(rshift(tmp0,i))){continue;}
					_mm_store_si128(idp,rshift(tmp0,i));
					Round_2();
				}
			}

			//初始化数组
			j0 = 0;
			memset(a0, 0, sizeof(a0));
			for (si8 i=0; i<N1; i++) { d0[i] = 1; f0[i] = i;} f0[N1] = N1;
			while (true)
			{	
				j0 = f0[0]; f0[0] = 0;
				if (j0==N1) break;
				else
				{
					idv[j0]=WtiForTravel[a0[j0]];
					tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));
					pr[N1]-=PDT_MaxProb[idv[j0]];
					a0[j0] = a0[j0] + d0[j0];
				}
				if ((a0[j0] == 0) || (a0[j0] == m0[j0]-1))
				{
					d0[j0] = -d0[j0];	f0[j0] = f0[j0+1];	f0[j0+1] = j0+1;
				}
				idv[j0]=WtiForTravel[a0[j0]];
				tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));
				pr[N1]+=PDT_MaxProb[idv[j0]];

				if(isValid(tmp0)&&((pr[N1]+Bn[Round-2])<(Bnc[Round-1]+1e-10))){
					r_pr[0]=pr[N1];
					_mm_store_si128(idp,tmp0);
					Round_2();
					for(int i=1;i<=7;i++){
						if(isDup(rshift(tmp0,i))){continue;}
						_mm_store_si128(idp,rshift(tmp0,i));
						Round_2();
					}
				}
			}



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

	fclose(fp_trails);
}