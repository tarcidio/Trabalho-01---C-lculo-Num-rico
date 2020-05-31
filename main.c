#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EULER 2.718281828459045235360287
#define MAX 100
#define precisao 0.0000000000000001

typedef double (*Funcao)(double);

double funcao_01(double x){
    return x - cos(x);
}

double derivada_primeira_funcao_01(double x){
    return 1 + sin(x);
}

double derivada_segunda_funcao_01(double x){
    return - cos(x);
}

double funcao_02(double x){
    return x*x*x - 9*x*x + 27*x - 27;
}

double derivada_primeira_funcao_02(double x){
    return 3*x*x - 18*x + 27;
}

double derivada_segunda_funcao_02(double x){
    return 6*x -18;
}

double funcao_03(double x){
    return pow(EULER,x) - cos(x);
}

double derivada_primeira_funcao_03(double x){
    return pow(EULER,x) + sin(x);
}

double derivada_segunda_funcao_03(double x){
    return pow(EULER,x) - cos(x);
}


double *metodo_bissecao(Funcao f, double esquerda_intervalo, double direita_intervalo, int *tamanho){
    double aproximacao;
    double novo;

    double* lista = (double*)malloc(1*sizeof(double));
    lista[0] = (esquerda_intervalo + direita_intervalo)/2;

    *tamanho = 1;
    if( (direita_intervalo - esquerda_intervalo) < precisao ){
        aproximacao = (direita_intervalo + esquerda_intervalo)/2;
        printf("\nNao houve necessidade da parte iterativa. Foi encontrada a raiz %.16lf", aproximacao);
        *tamanho = 0;
        return lista;
    }else{
        while(*tamanho <= MAX){
            novo = (esquerda_intervalo + direita_intervalo)/2;

            if( f(novo)*f(esquerda_intervalo) > 0 )
                esquerda_intervalo = novo;
            else
                direita_intervalo = novo;

            lista = realloc(lista, (*tamanho + 1)*sizeof(double));
            lista[*tamanho] = (direita_intervalo + esquerda_intervalo)/2;

            if( (fabs(direita_intervalo - esquerda_intervalo))/5 < precisao ){//por algum motivo, apenas a diferenca estava causando erro, entao apliquei a divisao
                aproximacao = (direita_intervalo + esquerda_intervalo)/2;
                printf("\nApos %d iteracoes, foi encontrada a raiz %.16lf", *tamanho, aproximacao);
                return lista;
            }

            (*tamanho)++;
        }
    }
    return NULL;
}

double *metodo_Newton(Funcao f, Funcao df, double aprox_inicial, int *tamanho){
        *tamanho = 1;
        double antigo, novo;

        double *lista = (double*)malloc(sizeof(double));
        lista[0] = aprox_inicial;

        if(fabs( f(aprox_inicial) ) < precisao){
            printf("\nNao houve necessidade da parte iterativa. Foi encontrada a raiz %.16lf", aprox_inicial);
            return lista;
        }else{
            antigo = aprox_inicial;
            while(*tamanho <= MAX){
                novo = antigo - f(antigo)/df(antigo);

                lista = realloc(lista, (*tamanho + 1)*sizeof(double));
                lista[*tamanho] = novo;

                if(fabs(f(novo)) < precisao || fabs(novo - antigo) < precisao){
                    printf("\nApos %d iteracoes, foi encontrada a raiz %.16lf", *tamanho, novo);
                    return lista;
                }
                antigo = novo;

                (*tamanho)++;
            }
        }
        return NULL;
}

double* metodo_Halley(Funcao f, Funcao df, Funcao ddf, double aprox_inicial, int *tamanho){
        *tamanho = 1;
        double antigo, novo;

        double *lista = (double*)malloc(sizeof(double));
        lista[0] = aprox_inicial;

        if(fabs(f(aprox_inicial)) < precisao){
            printf("\nNao houve necessidade da parte iterativa. Foi encontrada a raiz %.16lf", aprox_inicial);
            return lista;
        }else{
            antigo = aprox_inicial;
            while(*tamanho <= MAX){
                novo = antigo - f(antigo)*df(antigo)/( pow(df(antigo),2) - 0.5*f(antigo)*ddf(antigo));

                lista = realloc(lista, (*tamanho + 1)*sizeof(double));
                lista[*tamanho] = novo;

                if(fabs(f(novo)) < precisao || fabs(novo - antigo) < precisao){
                    printf("\nApos %d iteracoes, foi encontrada a raiz %.16lf", *tamanho, novo);
                    return lista;
                }
                antigo = novo;

                (*tamanho)++;
            }
        }
        return NULL;
}

double* metodo_secante(Funcao f, double aprox_inicial_01, double aprox_inicial_02,int *tamanho){
    *tamanho = 1;
    double antigo_01, antigo_02, novo;

    double *lista = (double*)malloc(sizeof(double));

    if(fabs(f(aprox_inicial_01)) < precisao){
        printf("\nNao houve necessidade da parte iterativa. Foi encontrada a raiz %.16lf", aprox_inicial_01);
        lista[0] = aprox_inicial_01;
        return lista;
    }else if(fabs(f(aprox_inicial_02)) < precisao || fabs(aprox_inicial_01 - aprox_inicial_02) < precisao){
        printf("\nNao houve necessidade da parte iterativa. Foi encontrada a raiz %.16lf", aprox_inicial_02);
        lista[0] = aprox_inicial_02;
        return lista;
    }else{
        antigo_01 = aprox_inicial_01;
        antigo_02 = aprox_inicial_02;
        while(*tamanho <= MAX){
            novo = antigo_02 - (f(antigo_02)*(antigo_02 - antigo_01))/(f(antigo_02) - f(antigo_01));

            lista = realloc(lista, (*tamanho + 1)*sizeof(double));
            lista[*tamanho] = novo;

            if(fabs(f(novo)) < precisao || fabs(novo - antigo_02) < precisao){
                printf("\nApos %d iteracoes, foi encontrada a raiz %.16lf", *tamanho, novo);
                return lista;
            }
            antigo_01 = antigo_02;
            antigo_02 = novo;

            (*tamanho)++;
        }
    }
    return NULL;
}

int main(void)
{
    int i, tamanho;


/**============================= PRIMEIRA EQUACAO ============================================= */


    printf("Primeira Equacao: x - cos(x) = 0");

    printf("\n\nMETODO DA BISSECAO (1 EQ.)");
    double* lista_bissecao_funcao_01 = metodo_bissecao(funcao_01, 0.6, 0.8, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_bissecao_funcao_01[i]);

    printf("\n\nMETODO DE NEWTON (1 EQ.)");
    double* lista_Newton_funcao_01  = metodo_Newton(funcao_01, derivada_primeira_funcao_01, 0.6, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Newton_funcao_01[i]);

    printf("\n\nMETODO DE HALLEY (1 EQ.)");
    double* lista_Halley_funcao_01  = metodo_Halley(funcao_01, derivada_primeira_funcao_01, derivada_segunda_funcao_01, 0.6, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Halley_funcao_01[i]);

    printf("\n\nMETODO DA SECANTE (1 EQ.)");
    double* lista_secante_funcao_01 = metodo_secante(funcao_01, 0.6, 0.8, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_secante_funcao_01[i]);


/**============================= SEGUNDA EQUACAO ============================================= */


    printf("\n\nSegunda Equacao: x^3-9x^2+27x-27=0");

    printf("\n\nMETODO DA BISSECAO (2 EQ.)");
    double* lista_bissecao_funcao_02 = metodo_bissecao(funcao_02, 2.6, 3.2, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_bissecao_funcao_02[i]);

    printf("\n\nMETODO DE NEWTON (2 EQ.)");
    double* lista_Newton_funcao_02  = metodo_Newton(funcao_02, derivada_primeira_funcao_02, 2.6, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Newton_funcao_02[i]);

    printf("\n\nMETODO DE HALLEY (2 EQ.)");
    double* lista_Halley_funcao_02  = metodo_Halley(funcao_02, derivada_primeira_funcao_02, derivada_segunda_funcao_02, 2.6, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Halley_funcao_02[i]);

    printf("\n\nMETODO DA SECANTE (2 EQ.)");
    double* lista_secante_funcao_02 = metodo_secante(funcao_02, 2.6, 3.2, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_secante_funcao_02[i]);


/**============================= TERCEIRA EQUACAO ============================================= */


    printf("\n\nTerceira Equacao: e^x-cos(x)=0");

    printf("\n\nMETODO DA BISSECAO (3 EQ.)");
    double* lista_bissecao_funcao_03 = metodo_bissecao(funcao_03, -0.95, 0.05, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_bissecao_funcao_03[i]);

    printf("\n\nMETODO DE NEWTON (3 EQ.)");
    double* lista_Newton_funcao_03  = metodo_Newton(funcao_03, derivada_primeira_funcao_03, 0.05, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Newton_funcao_03[i]);

    printf("\n\nMETODO DE HALLEY (3 EQ.)");
    double* lista_Halley_funcao_03  = metodo_Halley(funcao_03, derivada_primeira_funcao_03, derivada_segunda_funcao_03, 0.05, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_Halley_funcao_03[i]);

    printf("\n\nMETODO DA SECANTE (3 EQ.)");
    double* lista_secante_funcao_03 = metodo_secante(funcao_03, -0.95, 0.05, &tamanho);
    for(i = 0; i <= tamanho; i++)
        printf("\n%.16lf", lista_secante_funcao_03[i]);


    free(lista_bissecao_funcao_01);
    free(lista_Newton_funcao_01);
    free(lista_Halley_funcao_01);
    free(lista_secante_funcao_01);

    free(lista_bissecao_funcao_02);
    free(lista_Newton_funcao_02);
    free(lista_Halley_funcao_02);
    free(lista_secante_funcao_02);

    free(lista_bissecao_funcao_03);
    free(lista_Newton_funcao_03);
    free(lista_Halley_funcao_03);
    free(lista_secante_funcao_03);

    return 0;
    system("pause");
}
