/*
 * sorter.h
 *
 *  Created on: Oct 21, 2011
 *      Author: marius
 */

#ifndef SORTER_H_
#define SORTER_H_

#include "utils.h"

class Ssorter {
public:
	virtual ~Ssorter() {
	}
	virtual void sort(uint length, uint *data, sword *key) {
	}
};

class Sorter {
public:
	virtual ~Sorter() {
	}
	virtual void sort(uint length, uint *data, word *key) {
	}
};

class SMergeSorter: public Ssorter {
private:
	static const uint insertSortThreshold = 8;
	const uint maxBuffer;
	sword *destK;
	uint *destS;

public:
	SMergeSorter(uint maxBuffer) :
			maxBuffer(maxBuffer) {
		destS = new uint[maxBuffer];
		destK = new sword[maxBuffer];
	}
	~SMergeSorter() {
		delete[] destK;
		delete[] destS;
	}

	void sort(uint length, uint *sa, sword *key) {
		if (length <= insertSortThreshold) {
			sinsertSort(sa, length, key);
			return;
		}

		mergeSort(length, sa, key);
	}

	void mergeSort(uint length, uint *sa, sword*key) {
		uint s = 0;
		for (; s + insertSortThreshold < length; s += insertSortThreshold)
			sinsertSort(sa + s, insertSortThreshold, key + s);
		if (length - s > 1)
			sinsertSort(sa + s, length - s, key + s);

		// iterative merge sort
		uint *dS = destS;
		sword *dK = destK;
		for (uint bSize = insertSortThreshold; bSize < length; bSize <<= 1) {
			int n = 0;
			for (uint s = 0; s < length; s += (bSize << 1)) {
				uint i = s;
				uint m = s + bSize;
				if (m > length)
					m = length;
				uint j = m;
				uint e = m + bSize;
				if (e > length)
					e = length;
				for (; i < m && j < e; ++n) {
					if (key[i] <= key[j]) {
						dK[n] = key[i];
						dS[n] = sa[i];
						++i;
					} else {
						dK[n] = key[j];
						dS[n] = sa[j];
						++j;
					}
				}
				for (; i < m; ++n, ++i) {
					dK[n] = key[i];
					dS[n] = sa[i];
				}
				for (; j < e; ++n, ++j) {
					dK[n] = key[j];
					dS[n] = sa[j];
				}
			}
			uint *tmp = sa;
			sa = dS;
			dS = tmp;
			sword *tmpk = key;
			key = dK;
			dK = tmpk;
		}
		if (sa == destS) {
			memcpy(dK, key, length * sizeof(sword));
			memcpy(dS, sa, length * sizeof(uint));
		}
	}
};

class MergeSorter: public Sorter {
private:
	static const uint insertSortThreshold = 8;
	const uint maxBuffer;
	word *destK;
	uint *destS;

public:
	MergeSorter(uint maxBuffer) :
			maxBuffer(maxBuffer) {
		destS = new uint[maxBuffer];
		destK = new word[maxBuffer];
	}
	~MergeSorter() {
		delete[] destK;
		delete[] destS;
	}

	void sort(uint length, uint *sa, word *key) {
		if (length <= insertSortThreshold) {
			insertSort(sa, length, key);
			return;
		}

		mergeSort(length, sa, key);
	}

	void mergeSort(uint length, uint *sa, word*key) {
		uint s = 0;
		for (; s + insertSortThreshold < length; s += insertSortThreshold)
			insertSort(sa + s, insertSortThreshold, key + s);
		if (length - s > 1)
			insertSort(sa + s, length - s, key + s);

		// iterative merge sort
		uint *dS = destS;
		word *dK = destK;
		for (uint bSize = insertSortThreshold; bSize < length; bSize <<= 1) {
			int n = 0;
			for (uint s = 0; s < length; s += (bSize << 1)) {
				uint i = s;
				uint m = s + bSize;
				if (m > length)
					m = length;
				uint j = m;
				uint e = m + bSize;
				if (e > length)
					e = length;
				for (; i < m && j < e; ++n) {
					if (key[i] <= key[j]) {
						dK[n] = key[i];
						dS[n] = sa[i];
						++i;
					} else {
						dK[n] = key[j];
						dS[n] = sa[j];
						++j;
					}
				}
				for (; i < m; ++n, ++i) {
					dK[n] = key[i];
					dS[n] = sa[i];
				}
				for (; j < e; ++n, ++j) {
					dK[n] = key[j];
					dS[n] = sa[j];
				}
			}
			uint *tmp = sa;
			sa = dS;
			dS = tmp;
			word *tmpk = key;
			key = dK;
			dK = tmpk;
		}
		if (sa == destS) {
			memcpy(dK, key, length * sizeof(word));
			memcpy(dS, sa, length * sizeof(uint));
		}
	}
};

class RadixNotInPlaceSorter: public Sorter {
private:
	static const int charBuckets = 256;
	static const uint mergeSortThreshold = 1024;
	const uint passes;
	const uint maxBuffer;
	word *t_key;
	uint *t_sa;
	uint **bucketSizes;
	uint **bucketStarts;
	Sorter *ss;

	void sortCharNotInPlace(uint length, uint *sa, word *wkey, uint *tsa,
			word *tkey, int remainingKey) {
		const bool oddRound = (passes-remainingKey) & 1;
		if (length < mergeSortThreshold) {
			ss->sort(length, sa, wkey);
			if (oddRound) {
				memcpy(tsa, sa, length * sizeof(uint));
				memcpy(tkey, wkey, length * sizeof(word));
			}
//				for (uint i = 0; i < length; ++i) {
//					tsa[i] = sa[i];
//					tkey[i] = wkey[i];
//				}
			return;
		}

		int o = getMSDOffset(passes, remainingKey--);

		uint *bucketSize = bucketSizes[remainingKey];
		uint *bucketStart = bucketStarts[remainingKey];
		uchar *key = ((uchar*) wkey) + o;
		count(key, length, passes, bucketSize);
		bucketStart[0] = 0;
		prefixSum(bucketSize, bucketStart + 1, charBuckets - 1);

		for (uint i = 0; i < length; ++i, key += passes) {
			uint b = bucketStart[*key]++;
			tsa[b] = sa[i];
			tkey[b] = wkey[i];
		}
//		memcpy(sa, tsa, length * sizeof(uint));
//		memcpy(wkey, tkey, length * sizeof(word));

		if (remainingKey > 0) {
			key = ((uchar*) tkey) + o;
			for (uint i = 0; i < length;) {
				uchar c = key[i * passes];
				uint size = bucketSize[c];
				if (size > 1)
					sortCharNotInPlace(size, tsa + i, tkey + i, sa + i,
							wkey + i, remainingKey);
				else if (!oddRound) {
					sa[i] = tsa[i];
					wkey[i] = tkey[i];
				}
				i += size;
			}
		} else if (!oddRound) {
			memcpy(sa, tsa, length * sizeof(uint));
			memcpy(wkey, tkey, length * sizeof(word));
//			for (uint i = 0; i < length; ++i) {
//				sa[i] = tsa[i];
//				wkey[i] = tkey[i];
		}
	}

public:
	RadixNotInPlaceSorter(uint passes, uint maxBuffer) :
			passes(passes), maxBuffer(maxBuffer) {
		t_sa = new uint[maxBuffer];
		t_key = new word[maxBuffer];
		bucketSizes = create2Darray(passes, charBuckets);
		bucketStarts = create2Darray(passes, charBuckets);
		ss = new MergeSorter(mergeSortThreshold);
	}
	~RadixNotInPlaceSorter() {
		delete[] t_sa;
		delete[] t_key;
		delete2Darray(bucketSizes, passes);
		delete2Darray(bucketStarts, passes);
		delete ss;
	}

	void sort(uint length, uint *data, word *key) {
		sortCharNotInPlace(length, data, key, t_sa, t_key, passes);
	}

};
class SRadixLSD: public Ssorter {
private:
	static const int charBuckets = 256;
	static const int shortBuckets = 256 * 256;
	static const uint mergeSortThreshold = 4000;
	static const uint doubleCharSortThreshold = 40000;
	const uint passes;
	const uint maxBuffer;
	uint bucketSize[shortBuckets];
	uint bucketStart[shortBuckets];
	sword *tKey;
	uint *tSa;
	Ssorter *ss;

public:
	SRadixLSD(uint passes, uint maxBuffer) :
			passes(passes), maxBuffer(maxBuffer) {
		tSa = new uint[maxBuffer];
		tKey = new sword[maxBuffer];
		ss = new SMergeSorter(mergeSortThreshold);
	}
	~SRadixLSD() {
		delete[] tSa;
		delete[] tKey;
		delete ss;
	}

	void sort(uint length, uint *data, sword *wKey) {
		if (length < mergeSortThreshold) {
			ss->sort(length, data, wKey);
			return;
		}

		uint *sa = data;
		uint *t_sa = tSa;
		sword *wkey = wKey;
		sword *t_key = tKey;

		if (length > doubleCharSortThreshold) {
			int p = passes >> 1;
			for (uint r = p; r; --r) {
				int o = getLSDOffset(p, r);
				ushort *k = ((ushort*) wkey) + o;
				memset(bucketSize, 0, shortBuckets * sizeof(uint));
				count(k, length, p, bucketSize);

				if (exists(length, bucketSize, shortBuckets))
					continue;

				bucketStart[0] = 0;
				prefixSum(bucketSize, bucketStart + 1, shortBuckets - 1);

				uint * t = t_sa;
				t_sa = sa;
				sa = t;

				sword * tk = t_key;
				t_key = wkey;
				wkey = tk;

				for (uint i = 0; i < length; ++i, k += p) {
					uint b = bucketStart[*k]++;
					sa[b] = t_sa[i];
					wkey[b] = t_key[i];
				}
			}
		} else {
			for (uint r = passes; r; --r) {
				int o = getLSDOffset(passes, r);
				uchar *k = ((uchar*) wkey) + o;
				count(k, length, passes, bucketSize);

				if (exists(length, bucketSize, charBuckets))
					continue;

				bucketStart[0] = 0;
				prefixSum(bucketSize, bucketStart + 1, charBuckets - 1);

				uint * t = t_sa;
				t_sa = sa;
				sa = t;

				sword * tk = t_key;
				t_key = wkey;
				wkey = tk;

				for (uint i = 0; i < length; ++i, k += passes) {
					uint b = bucketStart[*k]++;
					sa[b] = t_sa[i];
					wkey[b] = t_key[i];
				}
			}
		}
		if (t_sa == data) {
			for (uint i = 0; i < length; ++i) {
				t_sa[i] = sa[i];
				t_key[i] = wkey[i];
			}
		}
	}

}
;
class RadixLSD: public Sorter {
private:
	static const int charBuckets = 256;
	static const int shortBuckets = 256 * 256;
	static const uint mergeSortThreshold = 1000;
	static const uint doubleCharSortThreshold = 40000;
	const uint passes;
	const uint maxBuffer;
	uint bucketSize[shortBuckets];
	uint bucketStart[shortBuckets];
	word *tKey;
	uint *tSa;
	Sorter *ss;

public:
	RadixLSD(uint passes, uint maxBuffer) :
			passes(passes), maxBuffer(maxBuffer) {
		tSa = new uint[maxBuffer];
		tKey = new word[maxBuffer];
		ss = new MergeSorter(mergeSortThreshold);
	}
	~RadixLSD() {
		delete[] tSa;
		delete[] tKey;
		delete ss;
	}

	void sort(uint length, uint *data, word *wKey) {
		if (length < mergeSortThreshold) {
			ss->sort(length, data, wKey);
			return;
		}

		uint *sa = data;
		uint *t_sa = tSa;
		word *wkey = wKey;
		word *t_key = tKey;

		if (length > doubleCharSortThreshold) {
			int p = passes >> 1;
			for (uint r = p; r; --r) {
				int o = getLSDOffset(p, r);
				ushort *k = ((ushort*) wkey) + o;
				memset(bucketSize, 0, shortBuckets * sizeof(uint));
				count(k, length, p, bucketSize);

				if (exists(length, bucketSize, shortBuckets))
					continue;

				bucketStart[0] = 0;
				prefixSum(bucketSize, bucketStart + 1, shortBuckets - 1);

				uint * t = t_sa;
				t_sa = sa;
				sa = t;

				word * tk = t_key;
				t_key = wkey;
				wkey = tk;

				for (uint i = 0; i < length; ++i, k += p) {
					uint b = bucketStart[*k]++;
					sa[b] = t_sa[i];
					wkey[b] = t_key[i];
				}
			}
		} else {
			for (uint r = passes; r; --r) {
				int o = getLSDOffset(passes, r);
				uchar *k = ((uchar*) wkey) + o;
				count(k, length, passes, bucketSize);

				if (exists(length, bucketSize, charBuckets))
					continue;

				bucketStart[0] = 0;
				prefixSum(bucketSize, bucketStart + 1, charBuckets - 1);

				uint * t = t_sa;
				t_sa = sa;
				sa = t;

				word * tk = t_key;
				t_key = wkey;
				wkey = tk;

				for (uint i = 0; i < length; ++i, k += passes) {
					uint b = bucketStart[*k]++;
					sa[b] = t_sa[i];
					wkey[b] = t_key[i];
				}
			}
		}
		if (t_sa == data) {
			for (uint i = 0; i < length; ++i) {
				t_sa[i] = sa[i];
				t_key[i] = wkey[i];
			}
		}
	}

}
;

class RadixMixedSorter: public Sorter {
private:
	static const int charBuckets = 256;
	static const uint insertSortThreshold = 16;
	static const uint mergeSortThreshold = 2048;
	static const uint lsdThreshold = 10000;
	const uint passes;
	const uint maxBuffer;
	word *tKey;
	uint *tSa;
	uint **bucketSizes;
	uint **bucketStarts;

public:
	RadixMixedSorter(uint passes, uint maxBuffer) :
			passes(passes), maxBuffer(maxBuffer) {
		tSa = new uint[maxBuffer];
		tKey = new word[maxBuffer];
		bucketSizes = create2Darray(passes, charBuckets);
		bucketStarts = create2Darray(passes, charBuckets);
	}

	~RadixMixedSorter() {
		delete[] tSa;
		delete[] tKey;
		delete2Darray(bucketSizes, passes);
		delete2Darray(bucketStarts, passes);
	}

	void sort(uint length, uint *data, word *wKey) {
		if (length < insertSortThreshold) {
			insertSort(data, length, wKey);
			return;
		}
		if (length < mergeSortThreshold) {
			mergeSort(length, data, wKey);
			return;
		}
		if (length < lsdThreshold) {
			lsdSort(length, data, wKey, 0);
			return;
		}
		notInPlaceMsdSort(length, data, wKey, passes);
	}

	void notInPlaceMsdSort(uint length, uint* data, word* wKey,
			int remainingKey) {
		int o = getMSDOffset(passes, remainingKey--);

		uint *bucketSize = bucketSizes[remainingKey];
		uint *bucketStart = bucketStarts[remainingKey];
		uchar *key = ((uchar*) wKey) + o;
		count(key, length, passes, bucketSize);
		bucketStart[0] = 0;
		prefixSum(bucketSize, bucketStart + 1, charBuckets - 1);

		memcpy(tSa, data, length * sizeof(uint));
		memcpy(tKey, wKey, length * sizeof(word));
		uchar *ct_key = ((uchar *) tKey) + o;
		uchar *cp = ct_key;
		for (uint i = 0; i < length; ++i, cp += passes) {
			uint b = bucketStart[*cp]++;
			data[b] = tSa[i];
			wKey[b] = tKey[i];
		}

		if (remainingKey > 0) {
			for (uint i = 0; i < length;) {
				uchar c = key[i * passes];
				uint size = bucketSize[c];
				if (size > 1) {
					uint *d = data + i;
					word *k = wKey + i;
					if (length < insertSortThreshold)
						insertSort(d, size, k);
					else if (length < mergeSortThreshold)
						mergeSort(size, d, k);
					else if (length < lsdThreshold)
						lsdSort(size, d, k, passes - remainingKey);
					else
						notInPlaceMsdSort(size, d, k, remainingKey);
				}
				i += size;
			}
		}
	}

	void lsdSort(uint length, uint *data, word *wKey, int start) {
		uint *sa = data;
		uint *t_sa = tSa;
		word *wkey = wKey;
		word *t_key = tKey;

		uint *bucketSize = bucketSizes[0];
		uint *bucketStart = bucketStarts[0];
		for (int r = passes; r > start; --r) {
			int o = getLSDOffset(passes, r);
			uchar *k = ((uchar*) wkey) + o;
			count(k, length, passes, bucketSize);

			if (exists(length, bucketSize, charBuckets))
				continue;

			int i;
			for (i = 0; i < charBuckets; ++i)
				if (length == bucketSize[i])
					break;

			if (i < charBuckets)
				continue;

			bucketStart[0] = 0;
			prefixSum(bucketSize, bucketStart + 1, charBuckets - 1);

			uint * t = t_sa;
			t_sa = sa;
			sa = t;

			word * tk = t_key;
			t_key = wkey;
			wkey = tk;

			for (uint i = 0; i < length; ++i, k += passes) {
				uint b = bucketStart[*k]++;
				sa[b] = t_sa[i];
				wkey[b] = t_key[i];
			}
		}
		if (t_sa == data) {
			for (uint i = 0; i < length; ++i) {
				t_sa[i] = sa[i];
				t_key[i] = wkey[i];
			}
		}
	}

	void mergeSort(uint length, uint *sa, word*key) {
		uint s = 0;
		for (; s + insertSortThreshold < length; s += insertSortThreshold)
			insertSort(sa + s, insertSortThreshold, key + s);
		if (length - s > 1)
			insertSort(sa + s, length - s, key + s);

		// iterative merge sort
		uint *dS = tSa;
		word *dK = tKey;
		for (uint bSize = insertSortThreshold; bSize < length; bSize <<= 1) {
			int n = 0;
			for (uint s = 0; s < length; s += (bSize << 1)) {
				uint i = s;
				uint m = s + bSize;
				if (m > length)
					m = length;
				uint j = m;
				uint e = m + bSize;
				if (e > length)
					e = length;
				for (; i < m && j < e; ++n) {
					if (key[i] <= key[j]) {
						dK[n] = key[i];
						dS[n] = sa[i];
						++i;
					} else {
						dK[n] = key[j];
						dS[n] = sa[j];
						++j;
					}
				}
				for (; i < m; ++n, ++i) {
					dK[n] = key[i];
					dS[n] = sa[i];
				}
				for (; j < e; ++n, ++j) {
					dK[n] = key[j];
					dS[n] = sa[j];
				}
			}
			uint *tmp = sa;
			sa = dS;
			dS = tmp;
			word *tmpk = key;
			key = dK;
			dK = tmpk;
		}
		if (sa == tSa) {
			memcpy(dK, key, length * sizeof(word));
			memcpy(dS, sa, length * sizeof(uint));
		}
	}
};
#endif /* SORTER_H_ */
