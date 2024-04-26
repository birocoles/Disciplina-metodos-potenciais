#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **aloca (FILE *arq, int grau);

double **libera (int linha, double **m);

void schmidt_norm_pnm (int grau, double latitude, double **pnm);

void geomagnetic_components(double longitude, double latitude, int grau, double razao_raio, double **pnm, double **gnm, double **hnm, double *X, double *Y, double *Z);

#define u0 1.25663706143592E-06

void main () {

	int i, j, k, l, grau, ordem;
	int n_latitude, n_longitude;
	double latitude, latitude_max, latitude_min, delta_latitude, v_min, delta_v, v;
	double longitude, longitude_min, longitude_max, delta_longitude, t_min, delta_t, t;
	double X, Y, Z;
	double raio_Terra, raio_calculo, razao_raio;
	double aux0, aux1;
	double **pnm, **derivada_pnm_colatitude, **gnm, **hnm;
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

		abort();

	}

	entrada = fopen(str, "r");

	fscanf (entrada, "%d", &grau);
	fscanf (entrada, "%lf %lf %d", &latitude_min, &latitude_max, &n_latitude);
	fscanf (entrada, "%lf %lf %d", &longitude_min, &longitude_max, &n_longitude);
	fscanf (entrada, "%lf %lf", &raio_Terra, &raio_calculo);
	
	ordem = grau;
	delta_latitude = (double)((latitude_max - latitude_min)/(n_latitude - 1));
	delta_longitude = (double)((longitude_max - longitude_min)/(n_longitude - 1));

	pnm = aloca (relatorio, grau);
	derivada_pnm_colatitude = aloca (relatorio, grau);
	gnm = aloca (relatorio, grau);
	hnm = aloca (relatorio, grau);

	fclose (entrada);

	/******** <== leitura do arquivo de input **************/

	/******** leitura do arquivo dos coeficientes de Gauss ==> **************/

	sprintf (str, "IGRF2005.txt", "r");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		abort();

	}

	entrada = fopen(str, "r");

	for (i = 1; i <= grau; i++) {
	
		for (j = 0; j <= i; j++) {
		
			fscanf (entrada, "%d %d %lf %lf", &k, &l, &aux0, &aux1);
			
			gnm[k][l] = aux0;
			hnm[k][l] = aux1;
		
		}
	
	}

	fclose (entrada);

	/******** <== leitura do arquivo dos coeficientes de Gauss **************/

	saida = fopen("output.txt", "w");
	
	fprintf (saida, " longitude   latitude                    X                    Y                    Z\n\n");

	//colatitude em radianos
   	v_min = (double)((latitude_min*3.14159265358979323846)/180.0);
   	v_min = (0.5*3.14159265358979323846) - v_min;
   	delta_v = (double)((delta_latitude*3.14159265358979323846)/180.0);

	//longitude em radianos
	t_min = (double)((longitude_min*3.14159265358979323846)/180.0);
   	delta_t = (double)((delta_longitude*3.14159265358979323846)/180.0);

	razao_raio = (double)(raio_Terra/raio_calculo);

	for (latitude = latitude_min, v = v_min; latitude <= latitude_max; latitude += delta_latitude, v += delta_v) {

		/* Cálculo dos harmônicos esféricos parcialmente normalizados
		de acordo com a normalização de Schmidt e da derivada numérica
		em relação à colatitude ==> */
		schmidt_norm_pnm (grau, (v + 0.000001), pnm);

		for (i = 0; i <= grau; i++) {
		
			for (j = 0; j <= i; j++) {
			
				derivada_pnm_colatitude[i][j] = pnm[i][j];
			
			}
		
		}
		
		schmidt_norm_pnm (grau, v, pnm);

		for (i = 0; i <= grau; i++) {
		
			for (j = 0; j <= i; j++) {
			
				derivada_pnm_colatitude[i][j] -= pnm[i][j];
				derivada_pnm_colatitude[i][j] *= 1000000.0;
			
			}
		
		}
		/* <== Cálculo dos harmônicos esféricos parcialmente normalizados
		de acordo com a normalização de Schmidt e da derivada numérica
		em relação à colatitude */

		for (longitude = longitude_min, t = t_min; longitude <= longitude_max; longitude += delta_longitude, t += delta_t) {

			geomagnetic_components(t, v, grau, razao_raio, pnm, derivada_pnm_colatitude, gnm, hnm, &X, &Y, &Z);

			fprintf (saida, "%10.5lf %10.5lf %20.5lf %20.5lf %20.5lf\n", longitude, latitude, X, Y, Z);

		}

	}

	fclose(saida);

	//Teste de validação da normalização de Schmidt	
	//for (i = 1; i <= 4; i++) {
	//
	//	for (j = 0; j <= i; j++) {
	//	
	//		printf ("%3d %3d %15.5lf %15.5lf\n", i, j, pnm[i][j], derivada_pnm_colatitude[i][j]);
	//	
	//	}
	//
	//}

	printf ("\n\nPrograma finalizado com sucesso!\n\n");

	system ("PAUSE");

}

/* Função que aloca a memória necessária para armazenar
os harmônicos esféricos e os coeficientes de Gauss */
double **aloca (FILE *arq, int grau) {

    double **m;  /* ponteiro para a matriz */
    int i;

    /* aloca as linhas da matriz */

    m = (double **)calloc((grau+1), sizeof(double *));

    if (m == NULL) {

        fprintf (arq, "Memoria Insuficiente (linhas)!\n\n");

        fclose (arq);

        system ("PAUSE");

		abort();

    }

    /* aloca as colunas da matriz */

    for (i = 0; i <= grau; i++) {

        m[i] = (double *)calloc((i+1), sizeof(double));

        if (m[i] == NULL) {

            fprintf (arq, "Memoria Insuficiente (colunas)!\n\n");

            fclose (arq);

            system ("PAUSE");

			abort();

        }

    }

    return (m); /* retorna o ponteiro para a matriz */

}

/* Função que desaloca a memória necessária para armazenar
os harmônicos esféricos e os coeficientes de Gauss */
double **libera (int linha, double **m) {

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

/* Função que calcula os harmônicos esféricos normalizados segundo
a normalização de Schmidt */
void schmidt_norm_pnm (int grau, double colatitude, double **pnm) {

	int i, j;
	double cosseno, seno;
	double aux0, aux1, aux2, aux3, aux4, aux5;

	cosseno = cos(colatitude);
	seno = sin(colatitude);

	/* Relação de recorrência para o cálculo dos harmônicos
	esféricos plenamente normalizados ==> */

    	/* Harmônico p00 */
	pnm[0][0] = 1.0;

    	/* Harmônico p10 */
	pnm[1][0] = cosseno*pow(3, 0.5);

	/* Harmônico p11 */
	pnm[1][1] = seno*pow(3, 0.5);

	/* Harmônico p20 */
	pnm[2][0] = (3*cosseno*cosseno) - 1.0;
	pnm[2][0] *= pow(1.25, 0.5);

    	/* Harmônico p21 */
    	pnm[2][1] = pow(15.0, 0.5)*cosseno*seno;

    	/* Harmônico p22 */
    	pnm[2][2] = pow(3.75, 0.5)*seno*seno;

	for(i = 3; i <= grau; i++) {

		aux0 = pow(((2*i) + 1), 0.5);
		aux1 = (double)(aux0/pow((2*i), 0.5));

		/* Harmônicos setoriais */
		pnm[i][i] = aux1*seno*pnm[i-1][i-1];

		/* Harmônicos pnn-1 */
		pnm[i][i-1] = aux0*cosseno*pnm[i-1][i-1];

        	/* Harmônicos pnn-2, pnn-3, pnn-4, ..., pn0 */
		aux2 = pow(((2*i) - 1.0), 0.5);
		aux3 = pow(((2*i) - 3.0), 0.5);

		for (j = i-2; j >= 0; j--) {

			aux4 = pow(((i+j)*(i-j)), 0.5);
			aux5 = pow(((i+j-1)*(i-j-1)), 0.5);

			pnm[i][j] = aux2*cosseno*pnm[i-1][j];
			pnm[i][j] -= ((double)(aux5/aux3))*pnm[i-2][j];
			pnm[i][j] *= (double)(aux0/aux4);

		}

	}
	/* <== Relação de recorrência para o cálculo dos harmônicos
	esféricos plenamente normalizados */

	/* transformação dos harmônicos esféricos plenamente normalizados
	para os harmônicos esféricos parcialmente normalizados de acordo
	com a normalização de Schmidt ==> */
	for (i = 0; i <= grau; i++) {
	
		aux0 = ((2.0*i) + 1.0);
		aux0 = pow(aux0, 0.5);
		aux0 = (double)(1.0/aux0);
	
		for (j = 0; j <= i; j++) {

			pnm[i][j] *= aux0;

		}
	
	}
	/* <== transformação dos harmônicos esféricos plenamente normalizados
	para os harmônicos esféricos parcialmente normalizados de acordo
	com a normalização de Schmidt */

}

/* Função que calcula as componentes X, Y e Z do campo geomagnético */
void geomagnetic_components(double longitude, double colatitude, int grau, double razao_raio, double **pnm, double **derivada_pnm_colatitude, double **gnm, double **hnm, double *X, double *Y, double *Z) {

	int i, j;
	double cosseno, seno;
	double aux0, aux1, aux2, aux3, aux4, aux5;

	(*X) = 0.0;
	(*Y) = 0.0;
	(*Z) = 0.0;
	
	aux4 = (double)(1.0/sin(colatitude));

	for (i = 1; i <= grau; i++) {

		aux0 = pow(razao_raio, (i+2));
		aux1 = 0.0;
		aux2 = 0.0;
		aux3 = 0.0;		

		for (j = 0; j <= i; j++) {

			cosseno = cos(j*longitude);
			seno = sin(j*longitude);
			aux5 = ((gnm[i][j]*cosseno) + (hnm[i][j]*seno));

			aux1 -= (aux5*derivada_pnm_colatitude[i][j]);
			aux2 += (((gnm[i][j]*seno) - (hnm[i][j]*cosseno))*j*pnm[i][j]);
			aux3 += (aux5*pnm[i][j]);

		}

		aux1 *= aux0;
		aux2 *= aux0;
		aux3 *= ((i + 1)*aux0);

		(*X) -= aux1;		
		(*Y) += aux2;
		(*Z) -= aux3;

	}
	
	(*Y) *= aux4;

}
