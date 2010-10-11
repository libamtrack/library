/*
 ============================================================================
 Name        : AT_D_RDD_GY.c
 Author      : Christoph Kolb
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "AmTrack.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[600];

	long n = 0;
	double r_m[500];
	double E_MeV_u=0;
	long particle_no=0;
	long material_no=0;
	long rdd_model=0;
	long rn = 0;
	double rdd_parameter[500];
	long er_model=0;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "r_m_input")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			char* nToken;
			nToken = strtok(token, " ");
			do {
				r_m[n] = atof(nToken);
				n++;
			} while (nToken = strtok(NULL, " "));
		}
		if (strstr(Text, "rdd_parameter")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			char* nToken;
			nToken = strtok(token, " ");
			do {
				rdd_parameter[rn] = atof(nToken);
				rn++;
			} while (nToken = strtok(NULL, " "));
		}
		if (strstr(Text, "E_Mev_u")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_MeV_u = atof(token);
		}
		if (strstr(Text, "er_model")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			er_model = atol(token);
		}
		if (strstr(Text, "rdd_model")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			rdd_model = atol(token);
		}
		if (strstr(Text, "material_no")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			material_no = atol(token);
		}
		if (strstr(Text, "particle_no")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			particle_no = atol(token);
		}
	}

	double D_RDD_Gy[n];

	AT_D_RDD_Gy(n, r_m, E_MeV_u, particle_no, material_no, rdd_model,
			rdd_parameter, er_model, D_RDD_Gy);
	char str[] = { "D_RDD_Gy:" };
	int i;
	for (i = 0; i < n; i++) {
		char text[1024];
		sprintf(text, " %f", D_RDD_Gy[i]);
		strcat(str, text);
	}
	strcat(str, "\n");
	fputs(str, f);

	if (f) {
		fclose(f);
	}

	return EXIT_SUCCESS;
}
