#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **aloca_pnm (FILE *arq, int grau);

double **libera_pnm (int linha, double **m);

double fully_norm_pnm (int grau, double colatitude, double **pnm);

int main () {

	int i, j, grau;
	double colatitude, colatitude_max, colatitude_min, delta_colatitude, v_min, delta_v, v;
	double **pnm;
	char str[100];

	FILE *relatorio, *entrada, *saida;

	relatorio = fopen ("relatorio.txt", "w");

	/******** leitura do arquivo de input ==> **************/

	sprintf (str, "input.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	fscanf (entrada, "%d %lf %lf %lf", &grau, &colatitude_min, &colatitude_max, &delta_colatitude);

	pnm = aloca_pnm (relatorio, grau);

	fclose (entrada);

	/******** <== leitura do arquivo de input **************/

	saida = fopen("saida.txt", "w");

   	v_min = (double)((colatitude_min*3.14159265358979323846)/180.0);
   	delta_v = (double)((delta_colatitude*3.14159265358979323846)/180.0);

	for (colatitude = colatitude_min, v = v_min; colatitude <= colatitude_max; colatitude += delta_colatitude, v += delta_v) {

		fully_norm_pnm(grau, v, pnm);

		fprintf (saida, "%10.3lf ", colatitude);

		printf ("%10.3lf\n", colatitude);

		fprintf (saida, "%15.5E ", pnm[0][0]);

        fprintf (saida, "%15.5E ", pnm[1][0]);

        fprintf (saida, "%15.5E ", pnm[1][1]);

		for (i = 2; i <= grau; i++) {

            fprintf (saida, "%15.5E ", pnm[i][0]);

			for (j = 1; j <= i; j++) {

				fprintf (saida, "%15.5E ", pnm[i][j]);

			}

		}

        fprintf (saida, "\n");

	}

	fclose(saida);

	printf ("\n\nPrograma finalizado com sucesso!\n\n");

	system ("PAUSE");

	return 0;

}

double **aloca_pnm (FILE *arq, int grau) {

    double **m;  /* ponteiro para a matriz */
    int i;

    /* aloca as linhas da matriz */

    m = (double **)calloc((grau+1), sizeof(double *));

    if (m == NULL) {

        fprintf (arq, "Memoria Insuficiente (linhas)!\n\n");

        fclose (arq);

        system ("PAUSE");
        return 0;

        return (NULL);

    }

    /* aloca as colunas da matriz */

    for (i = 0; i <= grau; i++) {

        m[i] = (double *)calloc((i+1), sizeof(double));

        if (m[i] == NULL) {

            fprintf (arq, "Memoria Insuficiente (colunas)!\n\n");

            fclose (arq);

            system ("PAUSE");
            return 0;

            return (NULL);

        }

    }

    return (m); /* retorna o ponteiro para a matriz */

}

double **libera_pnm (int linha, double **m) {

    int  i;  /* variavel auxiliar */

    if (m == NULL) {

        return (NULL);

    }

    for (i = 0; i <= linha; i++) {

        free (m[i]); /* libera as linhas da matriz */

    }

    free (m); /* libera a matriz */

    return (NULL); /* retorna um ponteiro nulo */

}

double fully_norm_pnm (int grau, double colatitude, double **pnm) {

	int i, j;
	double cosseno, seno;
	double aux0, aux1, aux2, aux3, aux4, aux5;

    /* Harmônico p00 */
	pnm[0][0] = 1.0;

	cosseno = cos(colatitude);
	seno = sin(colatitude);

    /* Harmônico p10 */
	pnm[1][0] = cosseno*pow(3, 0.5);

	/* Harmônico p11 */
	pnm[1][1] = seno*pow(3, 0.5);

	/* Harmônico p20 */
	pnm[2][0] = (3*cosseno*cosseno) - 1.0;
	pnm[2][0] *= pow(1.25, 0.5);

    /* Harmônico p21 */
    pnm[2][1] = pow(30, 0.5)*cosseno*seno;

    /* Harmônico p22 */
    pnm[2][2] = pow(22.5, 0.5)*seno*seno;

	for(i = 3; i <= grau; i++) {

		aux0 = pow(((2*i) + 1), 0.5);
		aux1 = (double)(aux0/(2*i));

		/* Harmônicos setoriais */
		pnm[i][i] = aux1*seno*pnm[i-1][i-1];

		/* Harmônicos pnn-1 */
		pnm[i][i-1] = aux0*cosseno*pnm[i-1][i-1];

        /* Harmônicos pnn-2, pnn-3, pnn-4, ..., pn0 */
		aux2 = pow(((2*i) - 1), 0.5);
		aux3 = pow(((2*i) - 3), 0.5);

		for (j = i-2; j >= 0; j--) {

			aux4 = pow(((i+j)*(i-j)), 0.5);
			aux5 = pow(((i+j-1)*(i-j-1)), 0.5);

			pnm[i][j] = aux2*cosseno*pnm[i-1][j];
			pnm[i][j] -= ((double)(aux5/aux3))*pnm[i-2][j];
			pnm[i][j] *= (double)(aux0/aux4);

		}

	}

	return 0;

}
