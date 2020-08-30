//MAP3121 - Metodos Numericos e Aplicacoes - Exercicio Programa 2
//NUSP - Aluno
//
//8943365 - Larissa Kimie Takayama
//9348985 - Thomas Palmeira Ferraz



#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

const int M_PI = 3.14159265359;

const double h_min = 1e-20; //situacoes normais de teste 0.001
const double h_max = 0.01; // situacoes normais de teste 1
const float p =4.0;
const double c = 2.0;

int R = 1450;
const double C1 = 10e-9;
const double  C2 = 100e-9;
const double L = 18e-3;
const double  E = 1.17391304;
const double  e_max= 8.1818;
double  G_a= (-50.0/66.0)*1e-3;
double  G_b= (-9.0/22.0)*1e-3;
double  G_c= 4.591*1e-3;

double * funcao_teste1 (double t, double *x, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = (1+pow((x[i] - t), 2));
    }
    return resposta;
}

double * multiplica_vetor_por_escalar (double *v, double escalar, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = escalar*v[i];
    }
    return resposta;
}
double * soma_vetor (double *v1, double *v2, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = v1[i] + v2[i];
    }
    return resposta;
}

double* soma_vetor (double *v1, double *v2, double *v3, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = v1[i] + v2[i] + v3[i];
    }
    return resposta;
}

double * soma_vetor (double *v1, double *v2, double *v3, double *v4, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = v1[i]+v2[i]+v3[i]+v4[i];
    }
    return resposta;
}

double * soma_vetor (double *v1, double *v2, double *v3, double *v4, double *v5, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = v1[i]+v2[i]+v3[i]+v4[i]+v5[i];
    }
    return resposta;
}

double* soma_vetor (double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *resposta, int s){
    for (int i = 0; i < s; i++){
        resposta[i] = v1[i]+v2[i]+v3[i]+v4[i]+v5[i]+v6[i];
    }
    return resposta;
}

double maxima_diferenca(double *x1, double *x2, int s){
    double dif = 0.0;
    double dif_i;
    for (int i = 0; i < s; i++){
        dif_i = fabs(x1[i] - x2[i]);
        if (dif_i > dif){
            dif = dif_i;
        }
    }
    return dif;

}

void guarda_vetor (vector<double *> *x,  double *y, int s){
    double *z = new double [s];
    for(int i = 0; i< s; i++){
        z[i] = y[i];
    }
    x->push_back(z);
}

//Função F(t, x(t)) para teste 2
double* funcao_teste2 (double t, double*x, double* F, int m) {
    double A[4][4];
    double soma;
    int i, j;

    //Definindo A
    A[0][0] = -2.0;
    A[0][1] = -1.0;
    A[0][2] = -1.0;
    A[0][3] = -2.0;
    A[1][0] = 1.0;
    A[1][1] = -2.0;
    A[1][2] = 2.0;
    A[1][3] = -1.0;
    A[2][0] = -1.0;
    A[2][1] = -2.0;
    A[2][2] = -2.0;
    A[2][3] = -1.0;
    A[3][0] = 2.0;
    A[3][1] = -1.0;
    A[3][2] = 1.0;
    A[3][3] = -2.0;

    //Multiplicando AX
    for (i=0; i<4; i++){
        soma = 0.0;
        for (j=0; j<4; j++){
            soma = soma + (A[i][j] * x[j]);
        }
        F[i] = soma;
    }

    return F;

}

//Função F(t, x(t)) para teste 3
double* funcao_teste3(double t,double* x, double* F, int m){
   int i, j;
   double** A = new double*[m];
   for (i = 0; i< m; i++){
        A[i] = new double [m];
   }
   double soma;

   //Contruindo A
   for (i=0; i<m; i++){
        for(j=0; j<m; j++){
            if(i==j){
                A[i][j] = -2.0;
            }
            else {
                if( j==(i+1) | i==(j+1) ){
                    A[i][j] = 1.0;
                }
                else {
                    A[i][j] = 0.0;
                }
            }
        }
   }

   //Multiplicando AX
   for (i=0; i< m; i++){
        soma = 0.0;
        for (j=0; j<m; j++){
            soma = soma + (A[i][j] * x[j]);
        }
        F[i] = soma;
    }
    for (i = 0; i< m; i++){
        delete(A[i]);
    }
    delete (A);
    return F;
}


//Função erro em cada ponto para teste 1
double erro_teste1 (double x_estim, double t){
    return fabs(x_estim - (t+1.0/(1.0-t)));
}


//Função erro em cada ponto para teste 2
double* erro_teste2 (double* x_estim, double t, double* erro){

    erro[0] = fabs(x_estim[0] - ( (exp(-1.0*t)*sin(t)) + (exp(-3.0*t)*cos(3*t)) ));
    erro[1] = fabs(x_estim[1] - ( (exp(-1.0*t)*cos(t)) + (exp(-3.0*t)*sin(3*t)) ));
    erro[2] = fabs(x_estim[2] - ( -1.0*(exp(-1.0*t)*sin(t)) + (exp(-3.0*t)*cos(3*t)) ));
    erro[3] = fabs(x_estim[3] - ( -1.0*(exp(-1.0*t)*cos(t)) + (exp(-3.0*t)*sin(3*t)) ));

    return erro;
}


//Função erro em cada ponto para teste 3
double* erro_teste3 (double* x_estim, double t, double* erro, int m){
    double lamb1, lamb2, y;
    int i;

    for (i=0; i<m; i++){
        y = i/(m+1.0);
        lamb1 = 2.0*(1.0-cos(M_PI/(m+1.0)));
        lamb2 = 2.0*(1.0-cos(m*M_PI/(m+1.0)));
        erro[i] = fabs(x_estim[i] - ( (exp(-1.0*lamb1*t)*sin(M_PI*y)) + (exp(-1.0*lamb2*t)*sin(m*M_PI*y)) ));
    }

    return erro;
}

double funcao_g (double v){
    if (v <= -1*e_max) return (G_c*v + e_max*(G_c - G_b) + E*(G_b - G_a));
    else if (v > -1*e_max and v <= -1*E) return (G_b*v + (G_b - G_a)*E);
    else if (v > -1*E and v < E) return (G_a*v);
    else if (v >= E and v < e_max) return (G_b*v + (G_a - G_b)*E);
    else return (G_c*v + e_max*(G_b - G_c) + E*(G_a - G_c));
}

double* funcao_circuito (double t,double* x, double* F, int m){
// Vc1, Vc2, IL
    F[0] = (1.0/(R*C1))*(x[1] - x[0]) + (-1.0/C1)*funcao_g(x[0]);
    F[1] = (1.0/(R*C2))*(x[0] - x[1]) + (1.0/C2)*x[2];
    F[2] = (-1.0/L)*x[1];
    return F;
}

void printa_vetor (string x, double * v, int m){
    cout << x;
    for (int i = 0; i< m; i++){
        cout << v[i] << " ";
    }
    cout << endl;
}

void printa_vetor (string x, double * v, int m, fstream saida){
    saida << x;
    for (int i = 0; i< m; i++){
        saida << v[i] << " ";
    }
    saida << endl;
}

void RKF45 (double h, double ti, double tf, double eps, double *x0, vector<double *> *x, vector<double> *t, vector<double *> *x_barra, double* funcao (double t, double *x, double *resposta, int s), unsigned int m){
    double t_atual = ti;
    double tau;
    double *k1 = new double [m];
    double *k2 = new double [m];
    double *k3 = new double [m];
    double *k4 = new double [m];
    double *k5 = new double [m];
    double *k6 = new double [m];
    double *temporario = new double [m];
    double *temporario2 = new double [m];
    double *temporario3 = new double [m];
    double *temporario4 = new double [m];
    double *temporario5 = new double [m];
    double *temporario6 = new double [m];
    double *temporario7 = new double [m];

    double *x_atual = new double [m];
    double *x_barra_atual = new double [m];
    bool preciso = false;
    x->push_back(x0);
    x_barra->push_back(x0);
    t->push_back (ti);
    ofstream saida;
    saida.open ("Circuito.txt");
    saida.precision(16);
	for(int i=0; t_atual < tf; i++){
        while (!preciso){
            k1 = multiplica_vetor_por_escalar(funcao(t_atual, x->at(i), temporario, m),h, k1, m);
            k2 = multiplica_vetor_por_escalar(funcao(t_atual + h/4.0, soma_vetor(x->at(i), multiplica_vetor_por_escalar(k1, 1.0/4.0, temporario,m), temporario2, m),temporario3,m),h,k2,m);
            k3 = multiplica_vetor_por_escalar(funcao(t_atual + 3.0*(h/8.0), soma_vetor(x->at(i), multiplica_vetor_por_escalar(k1, (3.0/32.0), temporario,m), multiplica_vetor_por_escalar(k2, (9.0/32.0), temporario2,m), temporario3,m), temporario4,m),h,k3,m);
            k4 = multiplica_vetor_por_escalar(funcao(t_atual + (12.0/13.0)*h, soma_vetor(x->at(i), multiplica_vetor_por_escalar(k1, (1932.0/2197.0), temporario,m), multiplica_vetor_por_escalar (k2, -7200.0/2197.0, temporario2,m), multiplica_vetor_por_escalar (k3, (7296.0/2197.0), temporario3,m), temporario4,m), temporario5,m),h, k4,m);
            k5 = multiplica_vetor_por_escalar(funcao(t_atual + h, soma_vetor(x->at(i),multiplica_vetor_por_escalar(k1, (439.0/216.0), temporario,m), multiplica_vetor_por_escalar(k2, -8.0, temporario2,m), multiplica_vetor_por_escalar(k3, (3680.0/513.0), temporario3,m), multiplica_vetor_por_escalar(k4, -845.0/4104.0, temporario4,m), temporario5,m), temporario6,m),h,k5,m);
            k6 = multiplica_vetor_por_escalar(funcao(t_atual + h/2.0, soma_vetor(x->at(i), multiplica_vetor_por_escalar(k1, -8.0/27.0, temporario,m), multiplica_vetor_por_escalar(k2, 2.0, temporario2,m), multiplica_vetor_por_escalar(k3, - 3544.0/2565.0, temporario3,m), multiplica_vetor_por_escalar(k4, (1859.0/4104.0), temporario4,m),multiplica_vetor_por_escalar(k5, -11.0/40.0, temporario5,m), temporario6,m), temporario7,m), h,k6,m);
            x_barra_atual = soma_vetor(x->at(i), multiplica_vetor_por_escalar(k1, 16.0/135.0, temporario,m), multiplica_vetor_por_escalar(k3, 6656.0/12825.0, temporario3,m), multiplica_vetor_por_escalar(k4, 28561.0/56430.0, temporario4,m), multiplica_vetor_por_escalar(k5, -9.0/50.0, temporario5,m),multiplica_vetor_por_escalar(k6, 2.0/55.0, temporario2,m),x_barra_atual,m);
            x_atual = soma_vetor(x->at(i),multiplica_vetor_por_escalar(k1, 25.0/216.0, temporario,m),multiplica_vetor_por_escalar(k3, 1408.0/2565.0, temporario3,m), multiplica_vetor_por_escalar(k4, 2197.0/4104.0, temporario4,m), multiplica_vetor_por_escalar(k5, -1.0/5.0, temporario5,m),x_atual,m);
            tau = maxima_diferenca(x_barra_atual,x_atual,m)/h;
            if (tau < eps){
                preciso = true;
                t_atual += h;
            }
            cout << "tau: " << tau << endl;
            h = h*pow (((h*eps)/(c*maxima_diferenca(x_barra_atual, x_atual, m))), 1.0/p);
            if (h > h_max) {
                h = h_max;
            }
            if (h < h_min) {
                h = h_min;
            }
            if (h > tf - t_atual) {
                h = tf - t_atual;
            }
            cout << "h: " << h << endl;

        }
        preciso = false;
        guarda_vetor (x, x_atual, m);
        t->push_back (t_atual);

        guarda_vetor (x_barra, x_barra_atual, m);
        cout << "i =" << i << endl;
	}
	saida.close();
	delete (k1);
	delete (k2);
	delete (k3);
	delete (k4);
	delete (k5);
	delete (k6);
	delete (temporario);
	delete (temporario2);
	delete (temporario3);
	delete (temporario4);
	delete (temporario5);
	delete (x_atual);
	delete (x_barra_atual);
}


int main()
{
    vector<double *> *x = new vector<double *>;
    vector<double> *t = new vector <double>;
    vector<double *> *x_barra = new vector<double *>;

    // Teste 2

    /*unsigned int m = 4;
    double *x0 = new double [m];
    x0[0]= 1.0;
    x0[1] = 1.0;
    x0[2] = 1.0;
    x0[3] = -1.0;
    RKF45(0.01, 0.0, 2, 1e-5, x0 , x, t, x_barra, &funcao_teste2, m);*/

    // Teste 1

    /*int m = 1;
    double *x0 = new double [m];
    x0 [0] = -18.95;
    RKF45(0.1, 1.05, 3, 1e-5, x0 , x, t, x_barra, &funcao_teste1, m);*/

    // Teste 3

    /*unsigned int m = 7;
    double *x0 = new double [m];
    for (int i = 0; i< m; i++){
        double y = ((double)i)/(1.0+m);
        x0[i] = sin (M_PI*y) + sin (m*M_PI*y);
    }
    RKF45(0.1, 0.0, 2, 1e-5, x0 , x, t, x_barra, &funcao_teste3, m);*/

    // Circuito de Chua
    int m = 3;
    double *x0 = new double [m];
    x0[0] = -0.5;
    x0[1] = - 0.2;
    x0[2] = 0.0;
    RKF45(0.01, 0.0, 0.05, 1e-1, x0 , x, t, x_barra, &funcao_circuito, m);

    //desalocacao
    delete (t);
    int s = x->size();
    for (int i = 0; i< s; i++){
        delete (x->at(i));
        delete (x_barra->at(i));
    }
    delete (x);
    delete (x_barra);

    return 0;
}
