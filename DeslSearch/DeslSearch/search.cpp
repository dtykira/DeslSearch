#include "search.h"

#define N0 (SBOX_NUMBER)
#define Ci(i) (C0[i+1]-1)
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

bool isValid(__m128i tmp){
	__m128i mask;
	mask=_mm_set1_epi16(0x0003);
	
	__m128i tmp0,tmp1,tmp2,tmp3,tmpCmp,tmpRes,Result;

	//����һ��__m128i�ļĴ��������Ĵ���_A�е�8��16bit��������_Count������ͬ���߼����ƣ�
	//��λ���ֵΪ0,r0=srl(_A0, _Count), r1=srl(_A1, _Count), ... r7=srl(_A7, _Count), 
	//shifting in zeros
	tmp0=_mm_srli_epi16(tmp, 4);// ||**** �� 0000||

	//����һ��__m128i�ļĴ�����r=srl(_A, _Imm * 8),   _Imm must be an immediate,  
	//shifting in zeros
	//print128(tmp);
	//print128(tmp0);
	tmp1=_mm_srli_si128(tmp0, 2);//print128(tmp1);
	tmp2=_mm_slli_si128(tmp0, 14);//print128(tmp2);
	tmp3=_mm_xor_si128(tmp1,tmp2);//print128(tmp3);
	
	tmpCmp=_mm_xor_si128(tmp,tmp3);//print128(tmpCmp);
	tmpRes=_mm_and_si128(tmpCmp,mask);
	
	u64 *p = (u64 *) &tmpRes;
	return (p[0]==0&&p[1]==0);
	return 0;
}

__m128i transform(u16 idv,si8 si){
	ALIGNED_TYPE_(u16,8) t[8];
	memset(t,0,16);
	t[si]=idv;
	__m128i *tmp;
	tmp=(__m128i *)t;
	return *tmp;
}

void Round_2(){
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od[0][si]);
	}fprintf(fp_trails,"\t%f\n",r_pr[0]);
}

void Round_1(){
	fp_trails=fopen("trails.txt","w");
	__m128i *idp;
	idp=(__m128i *)(r_od[0]);//idpָ��r_od
	__m128i tmp0;//��������ֵľֲ�����
	__m128i tmp;

	

	//0��S��
	//memset(r_od[0],0,SBOX_NUMBER*sizeof(u16));
	*idp=_mm_setzero_si128();
	firstRoundActive=false;
	Round_2();

	//1��S��
	firstRoundActive=true;
	for(int si=0;si<SBOX_NUMBER;si++){
		if(si!=0){r_od[0][si-1]=0;}
		for(int p=0;p<PR_NUMBER;p++){
			r_pr[0]=Prob[p];
			if((r_pr[0]+Bn[Round-2])>=(Bnc[Round-1]+1e-10)){break;}
			
			for(int i=0;i<PDT_1_Non0Num[p];i++){
				r_od[0][si]=PDT_1_Non0Val[p][i];
				if(!isValid(*idp)){continue;}
				//r_an,r_ai,r_pr
				Round_2();
			}
		}
	}

	//��ԾS������������
	prType pr[SBOX_NUMBER+1]={0};//���������
	// ��һ�������Ϊ������С�䶯�Ĵ��������ЩS�л�Ծ��
	// ����Ҫ��ע���� A0[] �� C0[]���������Ǹ����ı���
	// ���� A0[] �� 0 �� 1 �����У��� i ���ǻ�Ծ���� A[i] Ϊ 1������Ϊ 0
	// ���� C0[] �� �������У��� i ����ԾS�е��±����ͨ����ѯ C0[i] �õ�
	si8 A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
	// ��һ�������Ϊ������С�䶯�Ĵ�����������һ�ֵ������ֵģ��ڸ���weight��������һ�������ֵ������֣�
	// ���������ϵ��㷨�������ɻ�ϻ�ϵͳ�ĸ���Gray�루�Ƕ�Ԫ�ģ����µ�һ�� a0[] ��ǰһ�� a0[] �Ļ�����ֻ�䶯һ��Ԫ��
	// ����Ҫ��ע���� m0[], a0[], j0 
	// ���У�
	//     -- m0[] ���������У���ÿ����ԾS���ڸ���weight��������һ�������ֵ������ֵĸ�������ʼ��֮�󽫲��ٱ䶯
	//     -- a0[] ���������У���ÿ����ԾS���ڸ���weight��������һ�������ֵ������ֵ��±꣬ÿ��ѭ�����ᱻ�Ķ����µ�һ�� a0[] ��ǰһ�� a0[] �Ļ�����ֻ�䶯һ��Ԫ��
	//     -- j0 ����������ÿ��ѭ���б䶯���Ǹ�Ԫ�ص��±�
	si8 m0[SBOX_NUMBER], a0[SBOX_NUMBER], d0[SBOX_NUMBER], f0[SBOX_NUMBER+1];
	si8 j0;
	j0 = 0;
	si8 N1;
	//��������ֵľֲ�����
	si8 idv[SBOX_NUMBER];

	for(int M0=2;M0<=2;M0++){//��ԾS������������.M0<=N0
		N1=M0;
		for (si8 i=0; i<=(N0-M0); i++) A0[i] = 0;	for (si8 i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for (si8 i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1; T0[1] = 0; F0[N0] = N0 - M0 + 1; I0 = N0 - M0; L0 = N0;
		//��ʼ����ǰ���S�в���Ծ�������M0��S�л�Ծ
		
		do
		{
			//int pr=1;
			//���ɵ�һ��������
			tmp0=_mm_setzero_si128();
			for(si8 i=0;i<M0;i++){
				m0[i]=SBOX_INPUTS_NUMBER-1;
				idv[i]=WtiForTravel[0];
				tmp0=_mm_xor_si128(tmp0,transform(idv[i],Ci(i)));
				pr[i+1]=pr[i]+PDT_MaxProb[idv[i]];
				//if((pr[i+1]+Bn[Round-2])>=(Bnc[Round-1]+1e-10)){break;}
			}
			fprint128(tmp0);
			fprintf(fp_trails,"\tp:%f",pr[M0]);
			if(isValid(tmp0)&&(pr[M0]+Bn[Round-2])<(Bnc[Round-1]+1e-10)){fprintf(fp_trails,"\tpass");
				r_pr[0]=pr[M0];
				_mm_store_si128(idp,tmp0);
				//Round_2();
			}
			fprintf(fp_trails,"\n");

			//��ʼ������
			j0 = 0;
			memset(a0, 0, sizeof(a0));
			for (si8 i=0; i<N1; i++) { d0[i] = 1; f0[i] = i;m0[i]=i+1;} f0[N1] = N1;
			while (true)
			{	
				j0 = f0[0]; f0[0] = 0;
				if (j0==N1) break;
				else
				{
					idv[j0]=WtiForTravel[a0[j0]];
					tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));
					pr[M0]-=PDT_MaxProb[idv[j0]];
					a0[j0] = a0[j0] + d0[j0];
				}
				if ((a0[j0] == 0) || (a0[j0] == m0[j0]-1))
				{
					d0[j0] = -d0[j0];	f0[j0] = f0[j0+1];	f0[j0+1] = j0+1;
				}
				idv[j0]=WtiForTravel[a0[j0]];
				tmp0=_mm_xor_si128(tmp0,transform(idv[j0],Ci(j0)));
				pr[M0]+=PDT_MaxProb[idv[j0]];

				fprint128(tmp0);
				fprintf(fp_trails,"\tp:%f",pr[M0]);
				if(isValid(tmp0)&&((pr[M0]+Bn[Round-2])<(Bnc[Round-1]+1e-10))){fprintf(fp_trails,"\tpass");
					r_pr[0]=pr[M0];
					_mm_store_si128(idp,tmp0);
					Round_2();
				}
				fprintf(fp_trails,"\n");
			}


			//������һ��S�л�Ծģʽ
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