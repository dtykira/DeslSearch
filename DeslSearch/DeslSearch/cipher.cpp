#include "cipher.h"

//P�û������صı�
int PTable[32] = { 15, 6, 19, 20, 28, 11, 27, 16,
0, 14, 22, 25, 4, 17, 30, 9,
1, 7, 23, 13, 31, 26, 2, 8,
18, 12, 29, 5, 21, 10, 3, 24
};

//��չ�û�E�е��õĺ�������input��indexλ�ÿ�ʼȡ6������
void pick6(u8* output, u32 input, int index){
	if (index > 26){
		int carry = 32 - index;
		*output = ((input << (index - 26)) ^ (input >> (58 - index))) & 0x3f;
	}
	else{
		*output = (input >> (26 - index)) & 0x3f;
	}
}

//��չ�û�E
void Expansion(u8* output, u32 input){
	for (int i = 0; i < 8; i++){
		pick6(output + i, input, ((4 * i - 1) % 32+32)%32);
	}
}

//P�û�
u32 Permutation(u32 x){
	u32 y = 0;
	u32 t = 0x80000000;
	for (int i = 0; i<32; i++, t = t >> 1){
		if (PTable[i] - i>0) y |= (x << (PTable[i] - i))&t;
		else y |= (x >> (i - PTable[i])) & t;
	}
	return y;
}

//8��4��������ת����32���ش�
u32 SboxOutput2word(u8* input){
	u32 output=0;
	for(int i=0;i<7;i++){
		output+=*(input+i);
		output<<=4;
	}
	output+=*(input+7);
	return output;
}

//������S���������8��4��������
//�����E���������8��6��������
//ʵ�����Բ�
void PE(u8* output, u8* input){
	u32 temp=SboxOutput2word(input);
	temp=Permutation(temp);
	Expansion(output,temp);
}
