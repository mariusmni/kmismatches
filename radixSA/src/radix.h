/*
 * radix.h
 *
 *  Created on: Oct 21, 2011
 *      Author: marius
 */

#ifndef RADIX_H_
#define RADIX_H_

#include <math.h>
#define LITTLE_ENDIAN_FLAG 1
#define WORD64 1
//#define DEBUG 1

#include "utils.h"
#include "sorter.h"

#define BUCKET_CONST_SIZE 1

#ifdef WORD64
#define PASSES 8
#else
#define PASSES 4
#endif

class Radix {
private:
	uchar *originalInput;
	sword *input;
	uint *sa;
	uchar *bucketFlag;
	sword* sbuffer;
	word* buffer;
	uint *bucket;
	uint length;
	Sorter *sorter;
	Ssorter *ssorter;

	int alphaSize;
	uint bitsPerChar[maxChar];
	uint minBitsPerChar;
	uint maxBitsPerChar;
	uint charCode[maxChar];
	sword leftAlignedCharCode[maxChar];
	sword shiftedCode[maxChar];

#define ForNonSingleSubBucket(start, n, code) { uint _subLen = 1; \
		for (int _i = n - 1; _i >= 0; --_i) {\
			if (isBitSet(bucketFlag, start + _i)) { \
				uint _subStart = start + _i; \
				if (_subLen > 1) { \
						code; \
				} \
				_subLen = 1; \
			} else \
				++_subLen; }}

	void forEachNonSingletonBucket(uint start, uint n, uint helper,
			void(Radix::*h)(uint, uint, uint)) {
		uint bLen = 1;
		for (int i = n - 1; i >= 0; --i)
			if (isBitSet(bucketFlag, start + i)) {
				uint bStart = start + i;
				if (bLen > 1)
					(this->*h)(bStart, bLen, helper);
				bLen = 1;
			} else
				++bLen;
	}

	void forEachNonSingletonBucket(uint start, uint n, uint helper1,
			uint helper2, void(Radix::*h)(uint, uint, uint, uint)) {
		uint bLen = 1;
		for (int i = n - 1; i >= 0; --i)
			if (isBitSet(bucketFlag, start + i)) {
				uint bStart = start + i;
				if (bLen > 1)
					(this->*h)(bStart, bLen, helper1, helper2);
				bLen = 1;
			} else
				++bLen;
	}

	sword *prepareInput() {
		alphaSize
				= renameAlphabet(originalInput, length, charCode, bitsPerChar);
#ifdef DEBUG
		printf("String length %d alphabet size %d\n", length, alphaSize);
#endif

		minBitsPerChar = bitsPerSWord;
		maxBitsPerChar = 0;
		for (int i = 0; i < maxChar; ++i) {
			leftAlignedCharCode[i] = ((sword) charCode[i]) << (bitsPerSWord
					- bitsPerChar[i]);
			if (bitsPerChar[i] != 0) {
				if (bitsPerChar[i] < minBitsPerChar)
					minBitsPerChar = bitsPerChar[i];
				if (bitsPerChar[i] > maxBitsPerChar)
					maxBitsPerChar = bitsPerChar[i];
			}
		}

		if (minBitsPerChar == 0) {
			printf("ERROR: each code should have at least one bit\n");
			exit(0);
		}

		sword *inputw = new sword[IWORDSIZE + 1 + (length * bitsFor(alphaSize)
				/ bitsPerSWord)]; // one more dummy word as sentinel
		sword *input = (sword *) inputw;
		int bitsPerToken = bitsPerSWord;
		int remainBits = bitsPerToken;
		sword l = 0;
		int nWords = 0;
		for (uint i = 0; i < length; ++i) {
			uchar c = originalInput[i];
			uint code = charCode[c];
			uint codeLen = bitsPerChar[c];
			remainBits -= codeLen;
			if (remainBits >= 0) {
				l = (l << codeLen) | code;
			} else {
				int bitsInNextWord = -remainBits;
				int bitsAtEnd = codeLen - bitsInNextWord;
				l = (l << bitsAtEnd) | (code >> bitsInNextWord);
				input[nWords++] = l;
				l = GET_LAST(code, bitsInNextWord);
				remainBits = bitsPerToken - bitsInNextWord;
			}
		}
		if (remainBits < bitsPerToken) {
			input[nWords++] = l << remainBits;
		}
		input[nWords] = 0;

		return inputw;
	}

	void traverseInput(uint bitsAtATime, sword *codeShiftedByBitsAtATime,
			uint *helper, void(Radix::*h)(uint *, uint, sword)) {

		sword w = 0;
		for (int i = length - 1; i >= 0; --i) {
			int c = originalInput[i];
			w |= codeShiftedByBitsAtATime[c];
			w >>= bitsPerChar[c];
			(this->*h)(helper, i, w);
		}
	}

	void counter(uint *s, uint pos, sword data) {
		s[data]++;
	}

	void leftToRightSpreader(uint *bucketStart, uint pos, sword data) {
		sa[bucketStart[data]++] = pos;
	}

	void rightToLeftSpreader(uint *bucketEnd, uint pos, sword data) {
		sa[--bucketEnd[data]] = pos;
	}

	sword *shiftCode(int bitsAtATime) {
		for (int i = 0; i < maxChar; ++i)
			shiftedCode[i] = ((sword) charCode[i]) << bitsAtATime;
		return shiftedCode;
	}

	void firstPass(int bitsPerFirstPass, const int n, uint *bStart, uint *bSize) {
		shiftCode(bitsPerFirstPass);

		for (int i = 0; i < n; ++i)
			bSize[i] = 0;

		traverseInput(bitsPerFirstPass, shiftedCode, bSize, &Radix::counter);

		bStart[0] = 0;
		prefixSum(bSize, bStart + 1, n - 1);

		traverseInput(bitsPerFirstPass, shiftedCode, bStart,
				&Radix::leftToRightSpreader);

		for (int i = 0; i < n; ++i)
			bStart[i] -= bSize[i];
	}

	sword getSWordAtBit(int p) {
		int w = DIV_BY_WORDSIZE(p);
		int b = MOD_WORDSIZE(p);
		return b ? stitch(input[w], input[w + 1], bitsPerSWord - b, b)
				: input[w];
	}

	sword getSWordAt(uint s, uint bitsToSkip, uint dummy) {
		return getSWordAtBit(s * maxBitsPerChar + bitsToSkip);
	}

	void copyWords(uint* lo, int bLen, word* buffer, int bitsToSkip) {
		for (; bLen--; ++lo, ++buffer) {
			int skip = bitsToSkip;
			int bit = *lo * maxBitsPerChar;
#ifdef LITTLE_ENDIAN_FLAG
			for (int j = IWORDSIZE - 1; j >= 0; --j)
#else
			for (int j = 0; j < IWORDSIZE; ++j)
#endif
			{
				buffer->e[j] = getSWordAtBit(bit + skip);
				skip += bitsPerSWord;
			}
		}
	}

	void copyOneWordPerPosition(uint* lo, int bLen, word* buffer,
			int bitsToSkip) {
		for (int i = 0; i < bLen; ++i, ++buffer)
			buffer->e[0] = getSWordAt(lo[i], bitsToSkip, 0);
	}

	void updateBucketStart(sword *key, int length, uchar *bucketStart,
			int bStart) {
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				setBit(bucketStart, bStart + i);
	}

	void updateBucketStart(sword *key, int length, uchar *bucketStart,
			int bStart, uchar flag) {
		bucketStart += bStart;
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				bucketStart[i] = flag;
	}
	void updateBucketStart(word *key, int length, uchar *bucketStart,
			int bStart) {
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				setBit(bucketStart, bStart + i);
	}

	void updateBucketStart(word *key, int length, uchar *bucketStart,
			int bStart, uchar flag) {
		bucketStart += bStart;
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				bucketStart[i] = flag;
	}

	int secondPassHandle(uint bStart, uint bLen, uint bitsToSkip) {
		uint * lo = sa + bStart;
		//		uint D = bitsToSkip / maxBitsPerChar;
		//		uint oddBits = bitsToSkip - D * maxBitsPerChar;
		uint bitsPerCounter = 0;
		//		if (NotHandledPeriods
		//				== detectAndTreatPeriods(lo, D, bLen, bitsPerSWord, oddBits, 0,
		//						&Radix::getSWordAt, bitsPerCounter))
		//		printTime("Bucket");
		copyWords(lo, bLen, buffer, bitsToSkip);
		//		printTime("Copy");
		sorter->sort(bLen, lo, buffer);
		//		printTime("Sort");
		updateBucketStart(buffer, bLen, bucketFlag, bStart);
		return ((int) bitsPerSWord * IWORDSIZE) - ((int) bitsPerCounter);
	}

	uchar *createBucketFlags(int bitVectorSize, uint *start, uint n) {
		uchar* f = new uchar[bitVectorSize];
		memset(f, 0, bitVectorSize);
		for (uint i = 0; i < n; ++i)
			setBit(f, start[i]);
		return f;
	}

	int secondPass(uint bitsToSkip) {
		uint minBits = IWORDSIZE * bitsPerSWord;
		ForNonSingleSubBucket(
				0,
				length,
				minBits = min(minBits, secondPassHandle(_subStart, _subLen, bitsToSkip)));
		return minBits;
	}

	sword getBucketSWord(uint s1, uint D, uint bitsPerBucket) {
		if (s1 >= length) {
			return 0;
		} else {

#ifdef WORD64
			uint s2 = s1 + D;
			s1 = bucket[s1];
			if (s2 >= length)
				return ((sword) s1) << bitsPerBucket;

			s2 = bucket[s2];
			return (((sword) s1) << bitsPerBucket) | s2;
#else
			return bucket[s1];
#endif
		}
	}

	void copyBucketNumbers(uint *sa, uint D, uint D2, sword *buf, int bLen,
			int bitsPerBucket) {
		for (int i = 0; i < bLen; ++i)
			*buf++ = getBucketSWord(sa[i] + D, D2, bitsPerBucket);
	}

	void assignBucketNumbers(int bStart, int bLen) {
		uint *lo = sa + bStart;
		int currentBucket = 0;
		int firstBucket = bStart + 1;
		for (int i = 0; i < bLen; ++i) {
			if (bucketFlag[bStart++])
				currentBucket = firstBucket + i;
			bucket[lo[i]] = currentBucket;
		}
	}

	enum answer {
		HandledPeriodsFullFit, HandledPeriodsNoFullFit, NotHandledPeriods
	};

	void detectPeriods(uint *lo, uint bLen, uint halfD, uint &periodLength,
			int &maxNbPeriods) {
		periodLength = 0;
		maxNbPeriods = 0;
		int currentPeriodSize = 1;
		for (uint i = 1; i < bLen; ++i) {
			uint diff = lo[i - 1] - lo[i];
			if (diff <= halfD) { // period
				currentPeriodSize++;
				if (periodLength && diff != periodLength) {
					printf("Horror: period %d and %d\n", diff, periodLength);
					printSuffix(lo[i], 2 * halfD);
					printSuffix(lo[i - 1], 2 * halfD);
				}
				periodLength = diff;
			} else {
				if (currentPeriodSize > maxNbPeriods)
					maxNbPeriods = currentPeriodSize;
				currentPeriodSize = 1;
			}
		}
		if (currentPeriodSize > maxNbPeriods)
			maxNbPeriods = currentPeriodSize;
	}

	void detectPeriods2(uint *lo, uint bLen, uint D, uint &periodLength,
			int &maxNbPeriods) {
		periodLength = 0;
		maxNbPeriods = 0;
		int currentPeriodSize = 1;
		for (uint i = 1; i < bLen; ++i) {
			uint diff = lo[i - 1] - lo[i];
			if (diff <= D) { // period
				if (periodLength && diff != periodLength) {
					periodLength = 0;
					maxNbPeriods = 0;
					return; // no good
				}
				periodLength = diff;
				currentPeriodSize++;
			} else {
				if (currentPeriodSize > maxNbPeriods)
					maxNbPeriods = currentPeriodSize;
				currentPeriodSize = 1;
			}
		}
		if (currentPeriodSize > maxNbPeriods)
			maxNbPeriods = currentPeriodSize;
	}

	answer treatPeriods(uint *lo, uint D, uint bLen, int wordSizeInBits,
			uint helper1, uint helper2, sword(Radix::*wp)(uint, uint, uint),
			uint &bitsPerCounter, int maxNbPeriods, int periodLength) {
		const int maxc = 2 * maxNbPeriods - 1;
		bitsPerCounter = bitsFor(maxc);
		int restShift = bitsPerCounter + wordSizeInBits - bitsPerSWord;
		if (restShift < 0)
			restShift = 0;
		const int cShift = bitsPerSWord - bitsPerCounter;
		const sword maxcShifted = ((sword) maxc) << cShift;
		const sword countIncrement = ((sword) 1) << cShift;

		for (uint i = 0; i < bLen;) {
			int l = lo[i];
			sword leaderEnd = (this->*wp)(l + D - periodLength, helper1,
					helper2);
			sword rest = (this->*wp)(l + D, helper1, helper2);
			if (leaderEnd == rest)
				return NotHandledPeriods; // can happen with the copy from next optimization

			int leaderTypeS = leaderEnd < rest;
			sword c = (leaderTypeS ? maxcShifted : 0) | (rest >> restShift);
			sbuffer[i] = c;
			if (leaderTypeS)
				for (++i; i < bLen && lo[i - 1] - lo[i] == periodLength; ++i)
					c -= countIncrement, sbuffer[i] = c;
			else
				for (++i; i < bLen && lo[i - 1] - lo[i] == periodLength; ++i)
					c += countIncrement, sbuffer[i] = c;

		}
		return restShift ? HandledPeriodsNoFullFit : HandledPeriodsFullFit;
	}

	answer detectAndTreatPeriods(uint *lo, uint D, uint bLen,
			int wordSizeInBits, uint helper1, uint helper2, sword(Radix::*wp)(
					uint, uint, uint), uint &bitsPerCounter) {
		uint periodLength;
		int maxNbPeriods;
		detectPeriods2(lo, bLen, D, periodLength, maxNbPeriods);

		if (maxNbPeriods * periodLength > 2 * D)
			return treatPeriods(lo, D, bLen, wordSizeInBits, helper1, helper2,
					wp, bitsPerCounter, maxNbPeriods, periodLength);
		return NotHandledPeriods;
	}

	uint getBucketLength(int bStart) {
		for (int j = bStart;;)
			if (isBitSet(bucketFlag, ++j))
				return j - bStart;
		return 0;
	}

	void printBucketStats() {
		int n = 32;
		int bb[n];
		for (int i = 0; i < n; ++i)
			bb[i] = 0;

		int b = 0;
		int mb = 0;
		int singletons = 0;
		int prevBStart = 0;
		for (int i = 0; i <= length; ++i) {
			if (isBitSet(bucketFlag, i)) {
				if (b > mb) {
					mb = b;
				}
				//				/*
				for (int k = 0x10000, j = 0; k > 0; k >>= 1, ++j)
					if (b >= k) {
						bb[j]++;
					}
				//				 */
				if (b == 1) {
					singletons++;
				}

				b = 1;
				prevBStart = i;
			} else {
				++b;
			}
		}
		printf("LARGEST bucket %d\n", mb);
		//		/*
		for (int i = 0x10000, j = 0; i > 0; i >>= 1, ++j)
			if (bb[j] > 0)
				printf("buckets at least %d : %d\n", i, bb[j]);
		//		 */
		printf("singletons %d\n", singletons);
	}

	uint largestBuck() {
		uint b = 0;
		uint mb = 0;
		for (uint i = 0; i <= length; ++i)
			if (isBitSet(bucketFlag, i)) {
				if (b > mb)
					mb = b;
				b = 1;
			} else
				++b;
		return mb;
	}

	void printSuffix(int s, int len) {
		printf(" s=%d len=%d [", s, len);
		printCharData(originalInput + s, len);
		printf("]\n");
	}
	void printSuffixTranslated(int s, int len) {
		printf("Translated s=%d len=%d: ", s, len);
		for (int i = 0; i < len; ++i)
			printf("%3d ", charCode[originalInput[s + i]]);
		printf("\n");
	}

	void printLargestBucket(uint expectedSorted) {
		uint lb = largestBuck();
		printf("Largest %d\n", lb);
		uint bLen = 0;
		int bStart = 0;
		if (expectedSorted > 1000) {
			expectedSorted = 1000;
		}
		for (uint i = 0; i <= length; ++i) {
			if (isBitSet(bucketFlag, i)) {
				if (bLen > 1 && lb == bLen) {
					uint *lo = sa + bStart;
					uint s = *lo;
					printf("Bucket of size %d starts with", bLen);
					printSuffix(s, expectedSorted);
					break;
				}
				bLen = 1;
				bStart = i;
			} else {
				++bLen;
			}
		}
	}

	bool checkIncreasingSuffixes(int bStart, int bLen, uint expectedSorted) {
		for (int i = bStart + 1; i < bStart + bLen; ++i) {
			int j;
			for (j = 0; j < expectedSorted; ++j) {
				if (sa[i] + j >= length || sa[i - 1] + j >= length) {
					break;
				} else {
					if (originalInput[sa[i] + j]
							!= originalInput[sa[i - 1] + j])
						break;
				}
			}
			if (j < expectedSorted && originalInput[sa[i] + j]
					< originalInput[sa[i - 1] + j]
			//				&& originalInput[mySA[i]] != originalInput[mySA[i] + 1]
			) {
				printf(
						"OHO suffix %d before suffix %d but they differ at %d\n",
						sa[i - 1], sa[i], j);
				printf("OHO suffix %d: ", sa[i - 1]);
				printSuffix(sa[i - 1], expectedSorted);
				printSuffixTranslated(sa[i - 1], expectedSorted);
				printf("OHO suffix %d: ", sa[i]);
				printSuffix(sa[i], expectedSorted);
				printSuffixTranslated(sa[i], expectedSorted);

				return false;
				break;
			}
		}
		//
		return true;
	}

	/* To avoid using a special character for the terminal $, we
	 * solve any suffixes of the form 00...0$ first
	 */
	void fixEndingZeros() {
		int p = 0;
		for (uint i = length - 1; i && charCode[originalInput[i]] == 0; --i) {
			if (sa[p] != i) {
				printf("Error ending zeros!\n");
				exit(0);
			}
			setBit(bucketFlag, ++p);
		}
		//		printf("Ending zeros %d\n", p);
	}

	void updateSingletons(uchar *notSingleton, int bStart, int bLen) {
		uint *lo = sa + bStart;
		uchar *flo = bucketFlag + bStart;
		for (int j = 1; j <= bLen; ++j) {
			if (flo[j - 1] && flo[j])
				notSingleton[lo[j - 1]] = 0;
		}
	}

	bool less(uint s1, uint s2) {
		uint d = originalInput[s1];
		d -= originalInput[s2];
		if (d == 0)
			return bucket[s1 + 1] < bucket[s2 + 1];
		return d < 0;
	}

	answer handleBucketRepeatsOnly(uint bStart, uint bLen, uint D, uint D2,
			uint bitsPerBucket, uchar bucketStartFlag) {

		uint *lo = sa + bStart;
		uint bitsPerCounter;
		uint wordSizeInBits = bitsPerBucket;
#ifdef WORD64
		wordSizeInBits <<= 1;
#endif
		answer a = detectAndTreatPeriods(lo, D, bLen, wordSizeInBits, D2,
				bitsPerBucket, &Radix::getBucketSWord, bitsPerCounter);
		if (a != NotHandledPeriods) {
			ssorter->sort(bLen, lo, sbuffer);
			updateBucketStart(sbuffer, bLen, bucketFlag, bStart,
					bucketStartFlag);
			assignBucketNumbers(bStart, bLen);
		}
		return a;
	}

	answer handleBucketWithRepeats(uint bStart, uint bLen, uint D, uint D2,
			uint bitsPerBucket, uchar bucketStartFlag) {

		uint *lo = sa + bStart;
		uint bitsPerCounter;
		int wordSizeInBits = bitsPerBucket;
#ifdef WORD64
		wordSizeInBits <<= 1;
#endif
		answer a = detectAndTreatPeriods(lo, D, bLen, wordSizeInBits, D2,
				bitsPerBucket, &Radix::getBucketSWord, bitsPerCounter);
		if (a == NotHandledPeriods)
			copyBucketNumbers(lo, D, D2, sbuffer, bLen, bitsPerBucket);
		ssorter->sort(bLen, lo, sbuffer);
		updateBucketStart(sbuffer, bLen, bucketFlag, bStart, bucketStartFlag);
		assignBucketNumbers(bStart, bLen);
		return a;
	}
	bool adjacentSuffixes(uint p1, uint p2, uint n) {
		uint *a = sa + p1;
		uint *b = sa + p2;
		for (; n--;)
			if (*a++ != *b++ - 1)
				return false;
		return true;
	}

	void copyAdjacentSuffs(uint d, uint s, uint n) {
		uint *a = sa + d;
		uint *b = sa + s;
		for (; n--;)
			*a++ = *b++ - 1;
	}

	void copyFlags(uint d, uint s, uint n, uint f) {
		uchar *a = bucketFlag + d;
		uchar *b = bucketFlag + s;
		for (; n--; ++a)
			if (*b++)
				*a = f;
	}

	void finalTouches(uint expectedSorted) {
		uint nsingletons = 0;
		for (uint i = 0; i < length; ++i)
			if (bucketFlag[i] && bucketFlag[i + 1])
				++nsingletons;
		if (nsingletons == length)
			return;

		assignBucketNumbers(0, length);
		//		printTime("Bucket assignment        ");

		const bool useSingletons = nsingletons > (length >> 1);
		uchar *singleton = NULL;
		if (useSingletons) {
			singleton = new uchar[length];
			for (uint i = 0; i < length; ++i)
				singleton[i] = 1;
			for (uint i = 0; i < length; ++i)
				if (!bucketFlag[i] || !bucketFlag[i + 1])
					singleton[sa[i]] = 0;
		}
		//		printTime("Singletons");

		const int bitsPerBucket = bitsFor(length);
		const int allowedRepeats = 32;
		int totalAccess = length;
		do {
			uint prevBLen = 0;
			uint prevBStart = 0;
			uint nextBLen = 0;
			uint nextBStart = 0;
			uint bStart = 0;
			uint bLen = 0;
			bool copyFromNext = false;
			uint copyDepth = expectedSorted;
			int canTriple = 1;
			//				printBucketStats();				printLargestBucket(200);
			for (int i = length - 1; i >= 0; --i) {

				bool prevShouldCopy = false;
				prevBLen = 0;
				if (!useSingletons || !singleton[i]) {
					prevBStart = bucket[i] - 1;
					if (bucketFlag[prevBStart] <= allowedRepeats) {
						if (prevBLen == 0 || prevBStart != bStart) {
							prevBLen = getBucketLength(prevBStart);
							prevShouldCopy = prevBLen == bLen
									&& adjacentSuffixes(prevBStart, bStart,
											bLen);
						}
					}
				}

				if (bLen > 1) {
					totalAccess += bLen;
					int flag = bucketFlag[bStart];

					if (copyFromNext) {
						copyDepth++;
						//						printf("Copy depth is now %d\n", copyDepth);

						bool treatPeriods = copyDepth >= sa[bStart] - sa[bStart
								+ 1];
						copyAdjacentSuffs(bStart, nextBStart, bLen);
						copyFlags(bStart, nextBStart, bLen, flag + 1);
						assignBucketNumbers(bStart, bLen);

						if (treatPeriods) {
							ForNonSingleSubBucket(
									bStart,
									bLen,
									handleBucketRepeatsOnly(_subStart, _subLen, copyDepth, expectedSorted, bitsPerBucket, flag+1));
							copyDepth = 0;
						}
					} else {
						uint D = flag * expectedSorted;
						answer a = handleBucketWithRepeats(bStart, bLen, D,
								expectedSorted, bitsPerBucket, flag + 1);
						if (a == HandledPeriodsNoFullFit)
							canTriple = 0;
						copyDepth = (flag + 1) * expectedSorted;
					}
					//					if (useSingletons)
					//						updateSingletons(singleton, bStart, bLen);
					bucketFlag[bStart] = flag + 1;
				}

				nextBStart = bStart;
				nextBLen = bLen;
				if (bStart == prevBStart)
					bLen = getBucketLength(prevBStart);
				else
					bLen = prevBLen;
				bStart = prevBStart;
				copyFromNext = prevShouldCopy;
			}

//			cout << "PASS"<<endl;
			//			printTime("PASS");
			uint nBuckets = 0;
			for (uint i = 0; i < length; ++i) {
				uchar f = bucketFlag[i];
				if (f > 1)
					bucketFlag[i] = 1;
				if (f)
					nBuckets++;
			}
			if (nBuckets >= length)
				break;
			expectedSorted *= (2 + canTriple);

		} while (1);
//		cout << "Avg access per suffix " << (totalAccess * 1.0 / length)
//				<< endl;
		if (useSingletons)
			delete[] singleton;
	}

	void printHuTuckerCompression(const int n, uint *freq, uint *name,
			uint *codeLen, int bitsAtATime, int bitsPerChar) {
		int compressedSize = 0;
		for (int i = 0; i < n; ++i)
			compressedSize += freq[i] * codeLen[name[i]];

		int currentSize = length * bitsPerChar;
		printf(
				"Compressed length would be %d bits versus current %d ratio %.2f\n",
				compressedSize, currentSize, compressedSize * 1.0 / currentSize);
	}

	void sort(int n, uint *names, uint *key) {
		SRadixLSD s(PASSES, n);
		sword keyw[n];
		for (int i = 0; i < n; ++i)
			keyw[i] = key[i];
		s.sort(n, names, keyw);
		for (int i = 0; i < n; ++i)
			key[i] = keyw[i];

	}

	void printSortedChars(const int n, uint * freq, uint *chars) {
		int sum = 0;
		for (int i = 0; i < n; ++i)
			sum += freq[i];

		double P = 0;
		for (int i = 0; i < n; ++i) {
			double p = freq[i] * 1.0 / sum;
			printf("%3d. freq[%3d] = %d (%.2f%%)\n", i, chars[i], freq[i], p
					* 100);
			P += p * p;
		}
		printf("P = %.4lf, in theory we need %.4lf characters\n", P, 3 * log(
				length) / -log(P));
	}

	void printCharDistributionStats(int bitsAtATime, int bitsPerChar) {
		const int buckets = 1 << bitsAtATime;
		uint freq[buckets];
		uint chars[buckets];
		for (int i = 0; i < buckets; ++i)
			freq[i] = 0;
		traverseInput(bitsAtATime, shiftCode(bitsAtATime), freq,
				&Radix::counter);

		int n = collectNonZero(freq, freq, chars, buckets);

		uint code[n];
		uint codeLen[n];
		HuTucker(n, freq, code, codeLen);

		sort(n, chars, freq);

		printf("Chars:     ");
		printIntDataR(chars, n);

		printf("Freqs:     ");
		printIntDataR(freq, n);

		printf("Bit codes: ");
		printIntDataR(code, chars, n);

		printf("Code len:  ");
		printIntDataR(codeLen, chars, n);

		printHuTuckerCompression(n, freq, chars, codeLen, bitsAtATime,
				bitsPerChar);

		printSortedChars(n, freq, chars);
	}

	void checkSortedBucket(uint bStart, uint bLen, uint expectedSorted) {
		uint *lo = sa + bStart;
		for (int i = 1; i < bLen; ++i)
			if (memcmp(originalInput + lo[i - 1], originalInput + lo[i],
					expectedSorted) != 0) {
				printf("These should not be in the same bucket\n");
				printSuffix(lo[i - 1], expectedSorted);
				printSuffix(lo[i], expectedSorted);
			}

	}
public:
	Radix(uchar *input, uint n) :
		originalInput(input), length(n) {
	}

	uint* build() {
		resetTime();
		//		printTime("BEGIN");
		init();
		input = prepareInput();
		//		printTime("Prepared input");

		//		printCharDistributionStats(maxBitsPerChar, maxBitsPerChar);

		const int bitsPerFirstPass = (16 / maxBitsPerChar) * maxBitsPerChar;

		const int nBuckets = 1 << bitsPerFirstPass;
		uint bucketSize[nBuckets];
		uint bucketStart[nBuckets];
		sa = new uint[length];
		firstPass(bitsPerFirstPass, nBuckets, bucketStart, bucketSize);
		//		printTime("First passed");

		int bitVectorSize = length + 1; //(length + 100) >> 3;
		bucketFlag = createBucketFlags(bitVectorSize, bucketStart, nBuckets);
		setBit(bucketFlag, length);

		fixEndingZeros();
		//		printTime("Flags and ending zeros");

		uint maxBucket = max(bucketSize, nBuckets);
		buffer = new word[maxBucket];
		//		sorter = new RadixNotInPlaceSorter(IWORDSIZE * PASSES, maxBucket);
		sorter = new RadixLSD(IWORDSIZE * PASSES, maxBucket);
#ifdef DEBUG
		printf("MaxBucket %d\n", maxBucket);
#endif
		//		printBucketStats();
		/*
		 uint s[] {2619396, 2619396};
		 copyWords(s, 2, buffer, bitsPerFirstPass);
		 for (int i = 0; i < 2; ++i) {
		 printSuffix(s[i], 100);
		 printSuffixTranslated(s[i], 100);
		 printChars(buffer[i], maxBitsPerChar);
		 }
		 */
		int minBitsSortedInSecondPass = secondPass(bitsPerFirstPass);
		delete sorter;
		delete[] buffer;
		delete[] input;
		uint expectedSorted = (bitsPerFirstPass + minBitsSortedInSecondPass)
				/ maxBitsPerChar;
		//		ForNonSingleSubBucket(0, length,
		//				checkSortedBucket(_subStart, _subLen, expectedSorted));
		//		checkIncreasingSuffixes(0, length, expectedSorted);

		//		printTime("Further sorted");
#ifdef DEBUG
		printf("now sorted by at least %d chars\n", expectedSorted);
#endif

		//		exit(0);
		//		checkDepth(0, length);
		maxBucket = largestBuck();
#ifdef DEBUG
		printf("MaxBucket %d\n", maxBucket);
		//		printBucketStats();
		//		printLargestBucket(1000);
#endif
		bucket = new uint[length];
		sbuffer = new sword[maxBucket];
		//		ssorter = new RadixNotInPlaceSorter(PASSES, maxBucket);
		ssorter = new SRadixLSD(PASSES, maxBucket);
		finalTouches(expectedSorted);
		//		printTime("Final touches");

		//				checkIncreasingSuffixes(0, length, length);//100);
		delete[] sbuffer;
		delete[] bucket;
		delete[] bucketFlag;
		delete ssorter;
		return sa;
	}
};

#endif /* RADIX_H_ */
