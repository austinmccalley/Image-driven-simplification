#pragma once
#include <stdlib.h>

void sort(unsigned int* A, unsigned int* B, unsigned int* C, unsigned int sid, unsigned int eid) {
	unsigned int i;
	unsigned int* tempA, * tempB, * tempC;
	unsigned int current1, current2, current0;

	if (sid >= eid)
		return;
	sort(A, B, C, sid, (sid + eid) / 2);
	sort(A, B, C, (sid + eid) / 2 + 1, eid);
	tempA = (unsigned int*)malloc(sizeof(unsigned int) * (eid - sid + 1));
	tempB = (unsigned int*)malloc(sizeof(unsigned int) * (eid - sid + 1));
	tempC = (unsigned int*)malloc(sizeof(unsigned int) * (eid - sid + 1));
	for (i = 0; i < eid - sid + 1; i++) {
		tempA[i] = A[i + sid];
		tempB[i] = B[i + sid];
		tempC[i] = C[i + sid];
	}
	current1 = sid;
	current2 = (sid + eid) / 2 + 1;
	current0 = sid;
	while ((current1 <= (sid + eid) / 2) && (current2 <= eid)) {
		if (tempA[current1 - sid] < tempA[current2 - sid]) {
			A[current0] = tempA[current1 - sid];
			B[current0] = tempB[current1 - sid];
			C[current0] = tempC[current1 - sid];
			current1++;
		}
		else if (tempA[current1 - sid] > tempA[current2 - sid]) {
			A[current0] = tempA[current2 - sid];
			B[current0] = tempB[current2 - sid];
			C[current0] = tempC[current2 - sid];
			current2++;
		}
		else {
			if (tempB[current1 - sid] < tempB[current2 - sid]) {
				A[current0] = tempA[current1 - sid];
				B[current0] = tempB[current1 - sid];
				C[current0] = tempC[current1 - sid];
				current1++;
			}
			else {
				A[current0] = tempA[current2 - sid];
				B[current0] = tempB[current2 - sid];
				C[current0] = tempC[current2 - sid];
				current2++;
			}
		}
		current0++;
	}
	if (current1 <= (sid + eid) / 2) {
		for (i = current1; i <= (sid + eid) / 2; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}
	if (current2 <= eid) {
		for (i = current2; i <= eid; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0) {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0 / 3) {
		col[0] = 1.0;
		col[1] = 1.0 - percentage * 3.0;
		col[2] = 1.0 - percentage * 3.0;
	}
	else if (percentage <= 2.0 / 3) {
		col[0] = 1.0;
		col[1] = percentage * 3.0 - 1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0 / 3) {
		col[0] = 3.0 - percentage * 3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}