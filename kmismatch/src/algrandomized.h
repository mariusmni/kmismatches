#ifndef __algrandomized_h__
#define __algrandomized_h__ 1
#include "kutil.h"

void translate(int *a, int n, int *alpha, int r, int **code, int bitsPerChar,
		int *T) {
	for (int i = 0, j = 0; i < n; ++i) {
		int c = binarySearch(a[i], alpha, r);
		if (c < 0)
			c = r;
		for (int k = 0; k < bitsPerChar; ++k)
			T[j++] = code[c][k];
	}
}
/**
 * Counts for every position i in the text the number of mismatches between the
 * text and the pattern aligned at position i of the text. This is the randomized
 * version by Dr. Raj.
 *
 * Returns an array of the same size as the text.
 */
int *randomizedKMismatches(int *t, int n, int *p, int m, int k) {
	// alphabet and frequency
	int alpha[m];
	int freq[m];
	int r = alphaFreq(p, m, alpha, freq);

	const int bitsPerChar = 2 * (k + 1);
	int *code[r + 1];
	srand(123);
	for (int i = 0; i <= r; ++i) {
		code[i] = new int[bitsPerChar];
		memset(code[i], 0, bitsPerChar * sizeof(int));
	}

	//	for (int i = 0; i <= r; ++i)
	//		for (int j = 0; j < bitsPerChar; ++j)
	//			code[i][j] = rand() % 2;
	for (int i = 0; i <= r; ++i) {
		for (int j = 0; j < bitsPerChar; j += 2) {
			int p;
			while (code[i][p = rand() % bitsPerChar] == 1)
				;
			code[i][p] = 1;
		}
	}
	//	for (int i = 0; i <= r; ++i)
	//		print("Code ", code[i], bitsPerChar);

	int M = m * bitsPerChar;
	int *P = new int[M];
	translate(p, m, alpha, r, code, bitsPerChar, P);

	Convolution conv(3 * M - 1);
	conv.transform(P, M);

	int *T = new int[2 * M];
	int *C = new int[2 * M];
	int *matches = new int[n];
	memset(matches, 0, n * sizeof(int));
	for (int i = 0; i < n; i += m) {
		int size = min(n - i, 2 * m);
		translate(t + i, size, alpha, r, code, bitsPerChar, T);
		conv.convolveWithTransformed(T, bitsPerChar * size, M, C);
		for (int j = 0; j < size; ++j)
			matches[i + j] = C[j * bitsPerChar];

	}
	for (int i = 0; i <= r; ++i)
		delete[] code[i];

	for (int j = 0; j < n; ++j) {
		int s = matches[j];
		int eM = min(n - j, m);
		matches[j] = round(4.0 * s / bitsPerChar - eM);
		matches[j] = eM - matches[j]; //mismatches
	}
	delete []P;
	delete []T;
	delete []C;
	return matches;
}

/** Returns 1 if a matches b with at most k mismatches,
 * 0 if a matches b with more than k mismatches
 */
int countMismatches(int *a, int n, int *b, int m, int k) {
	while (n && m && k >= 0) {
		if (*a++ != *b++)
			--k;
		--n;
		--m;
	}
	return k >= 0 ? 1 : 0;
}

/** Returns an array of 0/1 values - 1 if the patern aligned at that position
 * matches with at most k mismatches, 0 if it matches with more than k mismatches.
 *
 * All positions for which the # of mismatches is within [k-epsilon, k+epsilon] are
 * confirmed in a brute force manner.
 */
int *randomizedKMismatches(int *t, int n, int *p, int m, int k, int epsilon) {
	int *m1 = randomizedKMismatches(t, n, p, m, k);
	int a = k - epsilon;
	int b = k + epsilon;
//	int byHand = 0;
	for (int i = 0; i < n; ++i) {
		if (m1[i] < a)
			m1[i] = 1;
		else if (m1[i] > b)
			m1[i] = 0;
		else {
			m1[i] = countMismatches(t + i, n - i, p, m, k);
//			byHand++;
		}
	}
//	cout << "RANDOMIZED checked manually " << byHand << endl;
	return m1;
}

#endif
