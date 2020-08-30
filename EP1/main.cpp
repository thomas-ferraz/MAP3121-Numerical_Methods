//MAP3121 - Metodos Numericos e Aplicacoes - Exercicio Programa 1
//NUSP - Aluno
//
//8943365 - Larissa Kimie Takayama
//9348985 - Thomas Palmeira Ferraz



#include <iostream>
#include <fstream>
#include <exception>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

typedef struct {
    double** A;
    int nlinhas, ncolunas;
} mat;
typedef struct {
    double* v;
    int n;
} vet;

mat le_matriz (){
    ifstream matriz_arquivo;
    string matriz_nome;
    cout << "Insira o titulo do arquivo da matriz ou zero para o padr�o:" << endl;
    cin >> matriz_nome;

    if (matriz_nome == "0"){
        matriz_nome = "76A_Completa_D_Matriz.txt";
    }
    
    cout << matriz_nome;
    
    try{
        matriz_arquivo.open (matriz_nome);
    }
    catch (exception& e)
    {
        cout << "Standard exception: " << e.what() << endl;
    }
    
    int nlinhas, ncolunas, nnulos;
    matriz_arquivo >> nlinhas >> ncolunas >> nnulos;
    double** matriz = new double*[ncolunas];

    for (int i = 0; i < ncolunas; i++)
        matriz[i] = new double [nlinhas];
    cout << nlinhas << " " << ncolunas << " " << nnulos << endl ;

    for (int i = 0; i < ncolunas; i++)
        for (int j = 0; j < nlinhas; j++)
            matriz [i][j] = 0.0;
    
    int a, b;
    double c;
    for (int k = 0; k < nnulos ; k++){
        matriz_arquivo >> a >> b >> c;
        cout << a << " " << b << " " << c << endl; // recebe linha coluna
        matriz[b][a] = c; // coloca [coluna][linha]
    }
    
    mat matriz_r;
    matriz_r.A = matriz;
    matriz_r.nlinhas = nlinhas;
    matriz_r.ncolunas = ncolunas;
    return matriz_r;
}

vet le_vetor (int nlinhas){
    ifstream arquivo;
    string nome;
    cout << "Insira o titulo do arquivo do vetor ou zero para o padr�o:" << endl;
    cin >> nome;
    
    if (nome == "0"){
        nome = "76A_Completa_D_VetorB.txt";
    }
    
    cout << nome;
    
    try{
        arquivo.open (nome);
    }
    catch (exception& e)
    {
        cout << "Standard exception: " << e.what() << endl;
    }
    
    double* vetor = new double[nlinhas];
    double x;
    for (int k = 0; k < nlinhas; k++){
        arquivo >> x;
        vetor [k] = x;
    }

    vet vetor_r;
    vetor_r.v = vetor;
    vetor_r.n = nlinhas;
    return vetor_r;
}

void delete_matriz (mat matriz){
    for (int i = 0; i < matriz.ncolunas; i++)
        delete (matriz.A[i]);
    delete (matriz.A);
}
void delete_vetor (vet vetor){
    delete (vetor.v);
}

//Printa um vetor
void imprime_saida_emtxt(double *v, int tamanho)
{
    ofstream saida;
    string nome = "resultado_teste76A.txt";
    saida.open (nome);
    saida.setf (ios::scientific);
    saida.precision(16);
    int i;
    for (i=0; i<tamanho; i++)
    {
        saida << v[i] << endl;
    }
    saida << endl;
    // EQM
    saida.close();
}

//Calcula a norma de um vetor
double norma(double* v, int m, int inicio) //recebe um vetor e seu tamanho
{
    int i;
    double soma=0, norma=0;
    for (i=inicio; i<m; i++)
    {
        soma = soma + (v[i]*v[i]); //soma quadrados
    }
    norma = sqrt(soma);
    return norma;
}

//Multiplicacao escalar por vetor
void multEscalar(double* v, int m, double escalar, double* resultado) //recebe vetor, tamanho dele e um escalar
{
    int i;
    for (i=0; i<m; i++)
    {
        resultado[i] = v[i] * escalar;
    }
}


//Produto escalar
double prodEscalar(double* v1, int m, double* v2, int inicio) //recebe dois vetores, seu tamanho e retorna o produto escalar
{
    int i;
    double produto=0;
    for(i=inicio; i<m; i++)
    {
        produto = produto + (v1[i] * v2[i]);
    }
    return produto;
}

//Subtracao de dois vetores
void subVetor(double* v1, int m, double* v2, double* subtracao, int inicio) //recebe vetor1, tamanho dos vetores, recebe vetor2 e coloca o resultado no vetor subtracao
{
    int i;
    for(i=inicio; i<m; i++)
    {
        subtracao[i]=v1[i]-v2[i];
    }
}


//Faz householder para uma matriz A e resultado b e retorna o R e b_barra
void houseHolder(double** A, double* b, int linhas, int colunas)
{
    int i, j, k;
    int delta=-1;
    double *w = new double [linhas];
    double *resParcial = new double [linhas];
    double resParcial2, resParcial3;
    for(i=0; i<colunas; i++)
    {
        if(A[i][i]>=0)
        {
            delta = 1;
        }

        for(k=0; k<i; k++)// otimiza��o, separar em dois for (0 a i) zera e de i a linhas copia.
        {
            w[k] = 0;//zera os elementos da coluna que est�o acima da diagonal principal
        }
        for (k=i; k < linhas; k++){
            w[k] = A[i][k]; //copia a coluna depois da sua diagonal principal
        }

        resParcial3 = delta*norma(w, linhas, i);
        w[i] = w[i] + resParcial3; //com essa mudanca o vetor a virou o vetor w

        for(j=i; j<colunas; j++) //aplicar o w em todas as colunas de A
        {
            resParcial2 = 0;
            if (norma(w, linhas, i)!=0)
            {
                resParcial2 = 2*prodEscalar(w, linhas, A[j], i)/prodEscalar(w, linhas, w, i);
            }
            multEscalar(w, linhas, resParcial2, resParcial);
            subVetor(A[j], linhas, resParcial, A[j], i);
        }
        resParcial2 = 0;
        if (norma(w, linhas, i) != 0)
        {
            resParcial2 = 2*prodEscalar(w, linhas, b, i)/prodEscalar(w, linhas, w, i);
        }
        multEscalar(w, linhas, resParcial2, resParcial);
        subVetor(b, linhas, resParcial, b, i);
    }
    delete (w);
    delete (resParcial);
}


//Resolve sistema linear "escalonado"
void resolveSistema(mat A, vet b, double* x)
{
    x[A.ncolunas-1] = b.v[A.ncolunas-1]/A.A[A.ncolunas-1][A.ncolunas-1];
    for (int i = A.ncolunas-2; i >= 0; i--){
        x [i] = 0.0;
        for (int j = i+1; j < A.ncolunas; j++){
            x[i] -= A.A[j][i]*x[j];
        }
        x[i] = (x[i]+b.v[i])/A.A[i][i];
    }
}

//Calcula fluxos de potencia
//Feito para rodar com B quadrada!
void fluxoPot(double** B, int tamanho, double* teta, double** P)
{
    //Achar os bjk
    int i, j;
    double soma;
    for(i=0; i<tamanho; i++)
    {
        soma = 0;
        for(j=0; j<tamanho; j++)
        {
            if(i!=j)
            {
                //aqui o P � construido de forma est�tica ent�o P[linha][coluna]
                //ao contrario de B que foi construido de forma dinamica e B[coluna][linha]
                P[j][i] = B[j][i]*(teta[i]-teta[j]);
                soma = soma + P[i][j];
            }
        }
        P[i][i] = soma;
    }

}

double erroQuadraticoMedio (double* b_barra, int linhas, int colunas){
    double soma = 0.0;
    for (int i=colunas; i < linhas; i++ ){
        soma += b_barra[i]*b_barra[i];
    }
    return sqrt(soma);
}

int reordena_linhas (mat matriz, vet v){
    int linhas = matriz.nlinhas;
    int colunas = matriz.ncolunas;
    double ** A = matriz.A;
    double * b = v.v;
    double ** A_ord = new double [colunas];
    for (int i=0; i < colunas; i++){
        A_ord[i] = new double [linhas];
    }
    int *indices_inicio = new int [linhas];
    int *indices_fim = new int [linhas];
    // percorre de forma nao continua na memoria
    for (int j=0; j < linhas; j++){
        int i;
        for (i=0; i < colunas or A[i][j] == 0.0; i++){
        }
        indices_inicio[j] = i;
        int i_fim = i;
        for (; i < colunas; i++){
            if (A[i][j] != 0.0){
                i_fim = i;
            }
        }
        indices_fim[j] = i_fim;
    }
    vector <pair <int,int>> inicios;
    for (int k=0; k< linhas; k++){
        inicios.push_back(make_pair (indices_inicio[k],k));
    }
    sort (inicios.begin(), inicios.end());


}

int main()
{
    mat A = le_matriz();
    vet b = le_vetor (A.nlinhas);
    cout << "vetor" << endl;

    houseHolder (A.A, b.v, A.nlinhas, A.ncolunas);
    cout << endl;
    cout << "householder" << endl;

    double *resp = new double [A.ncolunas];
    resolveSistema (A, b, resp);
    
    imprime_saida_emtxt(resp, A.ncolunas);
    // para pegar a coluna i da matriz basta pegar A.matriz [i]
    cout << "EQM = " << erroQuadraticoMedio(b.v, A.nlinhas, A.ncolunas) << endl;
    double ** Pot = new double* [A.ncolunas];

    //fluxoPot()
    delete_matriz (A);
    delete (resp);
    delete_vetor (b);

    return 0;
}