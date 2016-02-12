#ifndef __algabrahamson_h__
#define __algabrahamson_h__ 1
#include "kutil.h"

/**
 * Counts for every position i in the text the number of mismatches between the
 * text and the pattern aligned at position i of the text. This is the O(n sqrt(m log m))
 * method by Abrahamson which balances convolutions with brute force counting
 *
 * Returns an array of the same size as the text.
 */
int *abrahamsonCountMismatches(int *t, int n, int *p, int m) {
	// alphabet and frequency
	int alpha[m];
	int freq[m];
	int r = alphaFreq(p, m, alpha, freq);

	//sort characters by frequency in the pattern
	mypair a[r];
	for (int i = 0; i < r; ++i) {
		a[i].a = alpha[i];
		a[i].b = freq[i];
	}
	qsort(a, r, sizeof(mypair), pair_b_rcmp);

	// Take the most frequent A characters and count matches using convolutions
	int A = round(sqrt(m * 1.0 / (1+intLog(m))));
	int *matches = new int[n];
	memset(matches, 0, n * sizeof(int));
	int l = min(r, A);
	int chr[l];
	for (int i = 0; i < l; ++i)
		chr[i] = a[i].a;
	convolutions_fftw(chr, l, t, n, p, m, matches);

	// Take the remaining characters and count matches in a brute force manner
	if (A < r) {
		int offsets[r];
		memset(offsets, 0, r * sizeof(int));
		int *offset[r];
		for (int i = 0; i < r; ++i)
			offset[i] = new int[1+m / A];
		int remSize = r - A;
		int remAlphabet[remSize];
		for (int i = l; i < r; ++i)
			remAlphabet[i-l] = a[i].a;
		qsort(remAlphabet, remSize, sizeof(int), int_cmp);
		for (int i = 0; i < m; ++i) {
			int c = binarySearch(p[i], remAlphabet, remSize);
			if (c >= 0)
				offset[c][offsets[c]++] = i;
		}
		count(remAlphabet, remSize, t, n, offsets, offset, matches);

		for (int i = 0; i < r; ++i)
			delete[] offset[i];
	}
	for (int i = 0; i < n; ++i) {
		int k = min(m, n - i) - matches[i];
		matches[i] = k; // mismatches
	}
	return matches;
}

#endif
