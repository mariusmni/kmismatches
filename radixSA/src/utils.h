/*
 * utils.h
 *
 *  Created on: Oct 17, 2011
 *      Author: marius
 */

#ifndef MY_RADIX_UTILS_H_
#define MY_RADIX_UTILS_H_
#include <time.h>
#include <iostream>
using namespace std;

#define TOP(s, n, a)  (a = s[n-1])
#define POP(s, n, a)  (a = s[--n])
#define PUSH(s, n, a) (s[n++] = a)
#define TOP2(s, n, a, b)  {TOP(s,n-1,b); TOP(s,n-2,a);}
#define POP2(s, n, a, b)  {POP(s,n,b); POP(s,n,a);}
#define PUSH2(s, n, a, b) {PUSH(s,n,a); PUSH(s,n,b);}
#define TOP3(s, n, a, b, c)  {TOP(s, n, c); TOP(s,n-1,b); TOP(s,n-2,a);}
#define POP3(s, n, a, b, c)  {POP(s, n, c); POP(s,n,b); POP(s,n,a);}
#define PUSH3(s, n, a, b, c) {PUSH(s,n,a); PUSH(s,n,b); PUSH(s,n,c);}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define FORWARD 1
#define BACKWARD -1

typedef unsigned int uint;
typedef unsigned char uchar;
//typedef unsigned long long int ulong;
const int uintBits = 8 * sizeof(uint);

#ifdef WORD64
typedef unsigned long long sword;
#define DIV_BY_WORDSIZE(p) p >> 6
#define MOD_WORDSIZE(p) p & 63
#else
typedef uint sword;
#define DIV_BY_WORDSIZE(p) p >> 5
#define MOD_WORDSIZE(p) p & 31
#endif

#define IWORDSIZE 1

struct word {
	sword e[IWORDSIZE];
	inline bool operator<(const word& b) {
		int i;
		for (i = IWORDSIZE - 1; i > 0 && e[i] == b.e[i]; --i)
			;
		return e[i] < b.e[i];
	}
	inline bool operator<=(const word& b) {
		int i;
		for (i = IWORDSIZE - 1; i > 0 && e[i] == b.e[i]; --i)
			;
		return e[i] <= b.e[i];
	}
	inline bool operator==(const word& b) {
		for (int i = 0; i < IWORDSIZE; ++i)
			if (e[i] != b.e[i])
				return false;
		return true;
	}
	inline bool operator!=(const word& b) {
		for (int i = 0; i < IWORDSIZE; ++i)
			if (e[i] != b.e[i])
				return true;
		return false;
	}
};

const int bitsPerSWord = 8 * sizeof(sword);
static const int maxChar = 256;

#ifdef LITTLE_ENDIAN_FLAG
#define getBufferStart(buffer, len) ((buffer) + (len) - 1)
#define getMSDOffset(increment, remainingKey) ((remainingKey) - 1)
#define getLSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#else
#define getBufferStart(buffer, len) (buffer)
#define getMSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#define getLSDOffset(increment, remainingKey) ((remainingKey) - 1)
#endif

sword mask[bitsPerSWord][bitsPerSWord]; // mask of i 1's shifted to left j positions (e.g. mask[3][2] = 0..011100)
sword maskLeft[bitsPerSWord]; // mask of i 1's as most signif bits (e.g. mask[3] = 11100..0)
sword maskRight[bitsPerSWord]; // mask of i 1's  (e.g. mask[3] = 0...00111 = 7)
int bits[256];

void init() {
	for (int i = 0; i < bitsPerSWord; ++i) {
		sword ones = (((sword) 1) << i) - 1;
		maskRight[i] = ones;
		maskLeft[i] = ones << (bitsPerSWord - i);
		for (int j = 0; j < bitsPerSWord; ++j, ones <<= 1)
			mask[i][j] = ones;
	}

	bits[0] = 1;
	bits[1] = 1;
	for (int i = 2, b = 2; i < 256; ++b)
		for (int j = i * 2; i < j; ++i)
			bits[i] = b;
}

#define GET_LAST(w, b) ((w) & maskRight[b])
#define GET_FIRST(w, b) (((w) & maskLeft[b]) >> ((bitsPerSWord) - (b)))
#define GET_MIDDLE(w, b, o) (((w) & mask[b][o]) >> (o))

#define MAX(a, b) ((a) > (b) ? (a) : (b))

inline int max(int a, int b) {
	return a > b ? a : b;
}
inline int min(int a, int b) {
	return a < b ? a : b;
}

/**
 * Take the last1 bits from w1 and the first2 bits from w2 and form a new word
 */
inline sword stitch(sword w1, sword w2, int last1, int first2) {
	return (GET_LAST(w1, last1) << first2) | GET_FIRST(w2, first2);
}

void printChars(word w, uint bitsAtaTime) {
	int o = bitsPerSWord - bitsAtaTime;
	int c;
	do {
		c = GET_MIDDLE(w.e[1], bitsAtaTime, o);
		printf("%d ", c);
		o -= bitsAtaTime;
	} while (o >= 0);
	c = stitch(w.e[1], w.e[0], bitsAtaTime + o, -o);
	printf("%d ", c);
	o += bitsPerSWord;
	o -= bitsAtaTime;
	do {
		c = GET_MIDDLE(w.e[0], bitsAtaTime, o);
		printf("%d ", c);
		o -= bitsAtaTime;
	} while (o >= 0);
	printf("\n");
}

uchar unsetMask[] {~0x1, ~0x2, ~0x4, ~0x8, ~0x10, ~0x20, ~0x40, ~0x80};
uchar setMask[] {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};

void unSetBit(uchar *arr, int bitNumber) {
	int byte = bitNumber >> 3;
	int bit = bitNumber & 7;
	arr[byte] &= unsetMask[bit];
}

int bitsFor(int n) {
	int b = 0;
	if (n >= 0x10000) {
		b += 16;
		n >>= 16;
	}
	if (n >= 0x100) {
		b += 8;
		n >>= 8;
	}
	return b + bits[n];
}

/*
 int commonMSBits(word m, word n) {
 int b = 0;
 for (int i = 32; i > 0; i >>= 1)
 if ((m >> i) == (n >> i)) {
 b += i;
 } else {
 m >>= i;
 n >>= i;
 }
 return b + ((n & 1) == (m & 1));
 }*/

inline bool isBitSet(uchar *arr, int bitNumber) {
	return arr[bitNumber];
	//	int byte = bitNumber >> 3;
	//	int bit = bitNumber & 7;
	//	return arr[byte] & setMask[bit];
}

inline void setBit(uchar *arr, int bitNumber) {
	arr[bitNumber] = 1;
	//	int byte = bitNumber >> 3;
	//	int bit = bitNumber & 7;
	//	arr[byte] |= setMask[bit];
}

inline void swapInt(uint *a, uint *b) {
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

inline void swapSWords(sword *a, sword *b) {
	sword tmp = *a;
	*a = *b;
	*b = tmp;
}

inline void swapWords(word *a, word *b) {
	word tmp = *a;
	*a = *b;
	*b = tmp;
}

void swap(uchar *a, uchar *b, int size) {
	uchar tmp;
	for (; size--; ++a, ++b) {
		tmp = *a;
		*a = *b;
		*b = tmp;
	}
}
/*
 class Handle {
 public:
 virtual ~Handle() {
 }
 virtual void handle(uint pos, word data) {
 printf("NO NO baby\n");
 }
 };

 class Counter: public Handle {
 public:
 Counter(uint *size) {
 s = size;
 }
 void handle(uint pos, word data) {
 s[data]++;
 }
 private:
 uint *s;
 };

 class Spreader: public Handle {
 public:
 Spreader(uint *sa, uint *bucketEnd) {
 this->sa = sa;
 this->bucketEnd = bucketEnd;
 }
 void handle(uint pos, word data) {
 sa[--bucketEnd[data]] = pos;
 }
 private:
 uint *sa;
 uint *bucketEnd;
 };
 */
inline bool exists(uint e, uint *a, int n) {
	for (int i = 0; i < n; ++i)
		if (a[i] == e)
			return true;
	return false;
}

void prefixSum(uint* s, uint *d, int length) {
	d[0] = s[0];
	for (int i = 1; i < length; ++i)
		d[i] = d[i - 1] + s[i];
}

void prefixSum(uchar *s, uchar *d, int length) {
	*d++ = *s++;
	for (--length; length--; ++d, ++s)
		*d = *s + *(d - 1);
}

int renameAlphabet(uchar *in, uint length, uint *code, uint *codeLength) {
	const int maxAlpha = 256;
	for (int i = 0; i < maxAlpha; ++i)
		code[i] = 0;

	for (; length--;)
		code[*in++] = 1;

	code[0] -= 1;
	prefixSum(code, code, maxAlpha);

	int alphaSize = code[maxAlpha - 1] + 1;
	int bits = bitsFor(alphaSize);
	for (int i = 0; i < maxAlpha; ++i)
		codeLength[i] = bits;

	return alphaSize;
}

clock_t _startTime = clock();
clock_t _endTime;

void resetTime() {
	//#ifdef DEBUG
	_startTime = clock();
	//#endif
}

float getTime() {
	_endTime = clock();
	return (float) (_endTime - _startTime) / 1000; //CLOCKS_PER_SEC;
}
void printTime(char *msg) {
	//#ifdef DEBUG
	float seconds = getTime();
	cout << msg << " " << seconds << endl;
	resetTime();
	//#endif
}

int max(int *a, uint n) {
	int m = *a;
	for (; --n;)
		if (*++a > m)
			m = *a;
	return m;
}

uint max(uint *a, uint n) {
	uint m = *a;
	for (; --n;)
		if (*++a > m)
			m = *a;
	return m;
}

#ifdef WORD64
sword max(sword *a, uint n) {
	sword m = *a;
	for (; --n;)
	if (*++a > m)
	m = *a;
	return m;
}
#endif

/*
 word min(word *a, uint n) {
 word m = *a;
 for (; --n;)
 if (*++a < m)
 m = *a;
 return m;
 }*/

void count(uchar* key, int length, int increment, uint* bucketSize) {
	for (int i = 0; i < maxChar; ++i)
		bucketSize[i] = 0;
	for (; length--; key += increment)
		++bucketSize[*key];
}

void count(ushort* key, int length, int increment, uint* bucketSize) {
	for (; length--; key += increment)
		++bucketSize[*key];
}

uint **create2Darray(int rows, int cols) {
	uint **a = new uint*[rows];
	for (int i = 0; i < rows; ++i) {
		a[i] = new uint[cols];
	}
	return a;
}

void delete2Darray(uint **a, int rows) {
	for (int i = 0; i < rows; ++i) {
		delete[] a[i];
	}
	delete[] a;
}

void sinsertSort(uint* sa, uint length, sword* key) {
	if (length == 2) {
		if (key[1] < key[0]) {
			swapSWords(key, key + 1);
			swapInt(sa, sa + 1);
		}
		return;
	}

	for (uint i = 1; i < length; ++i) {
		sword kp = key[i];
		uint t = sa[i];
		uint j;
		for (j = i; j > 0 && kp < key[j - 1]; --j) {
			key[j] = key[j - 1];
			sa[j] = sa[j - 1];
		}
		key[j] = kp;
		sa[j] = t;
	}
}

void insertSort(uint* sa, uint length, word* key) {
	if (length == 2) {
		if (key[1] < key[0]) {
			swapWords(key, key + 1);
			swapInt(sa, sa + 1);
		}
		return;
	}

	for (uint i = 1; i < length; ++i) {
		word kp = key[i];
		uint t = sa[i];
		uint j;
		for (j = i; j > 0 && kp < key[j - 1]; --j) {
			key[j] = key[j - 1];
			sa[j] = sa[j - 1];
		}
		key[j] = kp;
		sa[j] = t;
	}
}

#define PRINT_LIMIT 500
void printCharData(uchar *array, int nKey) {
	int e = (nKey > PRINT_LIMIT) ? PRINT_LIMIT : nKey;
	for (int i = 0; i < e; ++i) {
		uchar c = array[i];
		if (c == '\n')
			c = 'N';
		if (c == '\r')
			c = 'n';
		if (c == '\t')
			c = 'T';
		if (c == 0)
			printf("0");
		else
			printf("%c", c);
	}
	if (e < nKey)
		printf("... %d more", nKey - e);
}

void printCharDataR(uchar *array, int nKey) {
	int s = (nKey > PRINT_LIMIT) ? nKey - PRINT_LIMIT : 0;
	for (int i = nKey - 1; i >= s; --i) {
		//		uchar c = array[i];
		//		if (c == '\n')
		//			c = 'N';
		//		if (c == '\r')
		//			c = 'n';
		//		if (c == 0)
		//			printf("  0(\0) ");
		//		else
		printf("%3d ", array[i]);
	}
	if (s > 0)
		printf("... %d more", s);
	printf("\n");
}

void printIntData(uint *array, int nKey) {
	int e = (nKey > PRINT_LIMIT) ? PRINT_LIMIT : nKey;
	for (int i = 0; i < e; ++i)
		printf("%d ", array[i]);
	if (e < nKey)
		printf("... %d more", nKey - e);
	printf("\n");
}

void printIntData(uint *array, uint *permutation, int nKey) {
	int e = (nKey > PRINT_LIMIT) ? PRINT_LIMIT : nKey;
	for (int i = 0; i < e; ++i)
		printf("%d ", array[permutation[i]]);
	if (e < nKey)
		printf("... %d more", nKey - e);
	printf("\n");
}

void printIntDataR(uint *array, int nKey) {
	int s = (nKey > PRINT_LIMIT) ? nKey - PRINT_LIMIT : 0;
	for (int i = nKey - 1; i >= s; --i)
		printf("%d ", array[i]);
	if (s > 0)
		printf("... %d more", s);
	printf("\n");
}

void printIntDataR(uint *array, uint *permutation, int nKey) {
	int s = (nKey > PRINT_LIMIT) ? nKey - PRINT_LIMIT : 0;
	for (int i = nKey - 1; i >= s; --i)
		printf("%d ", array[permutation[i]]);
	if (s > 0)
		printf("... %d more", s);
	printf("\n");
}

int collectNonZero(uint *src, uint *dest, uint *pos, uint len) {
	int n = 0;
	for (int i = 0; i < len; ++i)
		if (src[i] > 0)
			dest[n] = src[i], pos[n] = i, ++n;
	return n;
}

#define INFIN 0x7fffffff
void HuTuckerInitialTreeOld(int n, uint *f, int *parent) {
	int blocks = n;
	int twon = 2 * n;
	int i, j;
	int newGuy = blocks;
	bool *original = new bool[blocks];
	int *name = new int[blocks];
	int *freq = new int[blocks];
	for (i = 0; i < blocks; ++i)
		original[i] = true, freq[i] = f[i], name[i] = i;

	for (i = 0; i < twon; ++i)
		parent[i] = -1;

	int *lowFreqIndex = new int[blocks];
	int *lowFreq = new int[blocks];
	do {
		for (i = 0; i < blocks; ++i) {
			lowFreq[i] = INFIN;
			for (j = i - 1; j >= 0 && (j + 1 == i || !original[j + 1]); --j)
				if (freq[j] < lowFreq[i])
					lowFreq[i] = freq[j], lowFreqIndex[i] = j;
			for (j = i + 1; j < blocks && (j - 1 == i || !original[j - 1]); ++j)
				if (freq[j] < lowFreq[i])
					lowFreq[i] = freq[j], lowFreqIndex[i] = j;
		}

		for (i = 0; i < blocks; ++i) {
			for (j = i + 1; j < blocks; ++j)
				if (lowFreqIndex[i] == j && lowFreqIndex[j] == i)
					break;
			if (j < blocks)
				break;
		}

		if (i < blocks) {
			parent[name[i]] = newGuy;
			parent[name[j]] = newGuy;

			original[i] = false;
			freq[i] += freq[j];
			name[i] = newGuy;
			newGuy++;

			for (i = j + 1; i < blocks; ++i) {
				original[i - 1] = original[i];
				freq[i - 1] = freq[i];
				name[i - 1] = name[i];
			}
			blocks--;
		}
	} while (blocks > 1);
}

void HuTuckerInitialTree(const int n, uint *f, int *parent) {
	int twon = 2 * n;

	int i, j;
	int newGuy = n;

	bool original[n];
	int name[n];
	int freq[n];
	for (i = 0; i < n; ++i)
		original[i] = true, freq[i] = f[i], name[i] = i;

	for (i = 0; i < twon; ++i)
		parent[i] = -1;

	int lowPos[n];
	int low[n];
	int lowPosR[n];
	int lowR[n];
	int m = n;
	do {
		low[0] = INFIN;
		for (i = 1; i < m; ++i)
			if (original[i - 1] || freq[i - 1] < low[i - 1])
				low[i] = freq[i - 1], lowPos[i] = i - 1;
			else
				low[i] = low[i - 1], lowPos[i] = lowPos[i - 1];

		lowR[m - 1] = INFIN;
		for (i = m - 2; i >= 0; --i)
			if (original[i + 1] || freq[i + 1] < lowR[i + 1])
				lowR[i] = freq[i + 1], lowPosR[i] = i + 1;
			else
				lowR[i] = lowR[i + 1], lowPosR[i] = lowPosR[i + 1];

		for (i = 0; i < m; ++i)
			if (lowR[i] < low[i])
				low[i] = lowR[i], lowPos[i] = lowPosR[i];

		for (i = 0; i < m; ++i) {
			j = lowPos[i];
			if (i < j && i == lowPos[j]) {
				parent[name[i]] = newGuy;
				parent[name[j]] = newGuy;

				original[i] = false;
				freq[i] += freq[j];
				name[i] = newGuy;
				newGuy++;

				name[j] = -1;
			}
		}

		for (i = 0; name[i] == -1; ++i)
			;

		for (j = 0; i < m; ++i)
			if (name[i] != -1) {
				original[j] = original[i];
				freq[j] = freq[i];
				name[j] = name[i];
				++j;
			}

		m = j;

	} while (m > 1);
}

void HuTuckerDepth(int twon, int *parent, int *depth) {
	for (int i = 0; i < twon; ++i)
		depth[i] = 0;

	for (int i = twon - 1; i >= 0; --i)
		if (parent[i] != -1)
			depth[i] = depth[parent[i]] + 1;
}

int root(int node, int *parent) {
	while (parent[node] != -1)
		node = parent[node];
	return node;
}

void HuTuckerLastTree(int n, int *depth, int *parent) {
	const int twon = 2 * n;
	int maxDepth = max(depth, twon);

	for (int i = 0; i < twon; ++i)
		parent[i] = -1;

	int node[twon];
	int k = 0;
	int newGuy = n;
	do {
		if (maxDepth > 0)
			for (int i = 0; i < n; ++i)
				if (depth[i] == maxDepth) {
					int j;
					for (j = k; j > 0 && i < node[j - 1]; --j)
						node[j] = node[j - 1];
					node[j] = i;
					k++;
				}

		int j = 0;
		for (int i = 1; i < k; i += 2) {
			int r1 = root(node[i - 1], parent);
			int r2 = root(node[i], parent);
			parent[r1] = parent[r2] = newGuy++;
			node[j++] = node[i];
		}
		k = j;

		maxDepth--;
	} while (k > 1 || maxDepth > 0);
}

void HuTuckerMinLeaf(int n, int *parent, int *minLeaf) {
	int twon = 2 * n;
	int i;
	for (i = 0; i < n; ++i)
		minLeaf[i] = i;

	for (; i < twon; ++i)
		minLeaf[i] = INFIN;

	for (i = 0; i < twon; ++i)
		if (parent[i] != -1)
			if (minLeaf[i] < minLeaf[parent[i]])
				minLeaf[parent[i]] = minLeaf[i];
}

void HuTuckerComputeCodes(int n, int *parent, uint *code, uint *codeLen) {
	const int twon = 2 * n;

	int minLeaf[twon];
	HuTuckerMinLeaf(n, parent, minLeaf);

	for (int i = 0; i < twon; ++i)
		code[i] = codeLen[i] = 0;

	for (int i = twon - 1; i >= 0; --i)
		if (parent[i] != -1) {
			int rightChild = minLeaf[i] > minLeaf[parent[i]];
			code[i] = (code[parent[i]] << 1) | rightChild;
			codeLen[i] = codeLen[parent[i]] + 1;
		}
}

void HuTucker(int n, uint *f, uint *code, uint *codeLen) {
	int twon = 2 * n;

	int parent[twon];
	HuTuckerInitialTree(n, f, parent);

	int depth[twon];
	HuTuckerDepth(twon, parent, depth);

	HuTuckerLastTree(n, depth, parent);

	uint code2[twon];
	uint codeLen2[twon];
	HuTuckerComputeCodes(n, parent, code2, codeLen2);

	for (int i = 0; i < n; ++i)
		code[i] = code2[i], codeLen[i] = codeLen2[i];
}

int renameAlphabetHuTucker(uchar *in, uint length, uint *code, uint *codeLen) {
	const int maxAlpha = 256;
	uint freq[maxAlpha];
	count(in, length, 1, freq);

	uint name[maxAlpha];
	int n = collectNonZero(freq, freq, name, maxAlpha);

	uint encod[n];
	uint encodLen[n];
	HuTucker(n, freq, encod, encodLen);

	for (int i = 0; i < maxAlpha; ++i)
		code[i] = 0, codeLen[i] = 0;

	for (int i = 0; i < n; ++i)
		code[name[i]] = encod[i], codeLen[name[i]] = encodLen[i];
	return n;
}

int main4() {
	uint freq[] {1, 2, 23, 4, 3, 3, 5, 19};
	int n = sizeof(freq) / sizeof(uint);
	printf("Freqs:     ");
	printIntData(freq, n);

	uint *code = new uint[n];
	uint *codeLen = new uint[n];
	HuTucker(n, freq, code, codeLen);

	printf("Bit codes: ");
	printIntData(code, n);
	printf("Code len:  ");
	printIntData(codeLen, n);

}
#endif
