/*
 * runtimes.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: marius
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "radix.h"


uchar* readData(const char * const filename, uint& n) {
	struct stat fileInfo;
	FILE *file;
	if (stat(filename, &fileInfo)) {
		printf("Unable to get stat of file %s \n", filename);
		return NULL;
	}
	n = fileInfo.st_size;
	uchar *result = new uchar[n];
	if (!(file = fopen(filename, "r"))) {
		printf("Unable to open file %s \n", filename);
		return NULL;
	}
	rewind(file);
	if (n > fread(result, sizeof(uchar), n, file)) {
		printf("Error reading file %s \n", filename);
		fclose(file);
		delete[] result;
		return NULL;
	}
	fclose(file);
	return result;
}

void printIntData(uint *data, uint n, char *outputFile) {
	FILE *f = fopen(outputFile, "w");
	for (uint i = 0; i < n; ++i) {
		fprintf(f, "%d ", data[i]);
	}
	fprintf(f, "\n");
	fclose(f);
}

int main(int argc, char *argv[]) {
	int i =1;
	char *c = (char *)&i;
#ifdef LITTLE_ENDIAN_FLAG
	bool fail = (*c == 0);
#else
	bool fail = (*c == 1);
#endif
	if (fail) {
		printf("Failed to detect endianness! Please manually set/unset LITTLE_ENDIAN_FLAG and recompile\n");
		exit(0);
	}



	if (argc < 3) {
		printf("Expected 2 arguments: input_file output_file\n");
		return 0;
	}
	char *inputFile = argv[1];
	char *outputFile = argv[2];
	uint n = 0;
	uchar *ustr = readData(inputFile, n);
	if (ustr == NULL)
		return 0;

	if (ustr[n - 1] == 10) { // remove line feed
		ustr[--n] = 0;
	}

	printf("===============================================\n");
	printf("string length %d of file %s \n", n, inputFile);
	clock_t _startTime = clock();
	Radix *r = new Radix(ustr, n);
	unsigned int *radixSA = r->build();
	delete r;
	clock_t _endTime = clock();
	float seconds = (float) (_endTime - _startTime) / 1000; //CLOCKS_PER_SEC;
	printf("RadixSA took [%.2fms]\n", seconds);
	printIntData(radixSA, n, outputFile);
	delete[] radixSA;

	exit(0);
	return 0;
}
