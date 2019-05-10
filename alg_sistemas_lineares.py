
import numpy as np
import matplotlib.pyplot as plt
import time

#matriz[linha][coluna]


#Le as matrizes do arquivo txt
def lerArquivo(nomeArquivo):
    arq = open(nomeArquivo, 'r')
    texto = arq.readlines()
    matriz = []

    for linha in texto:
        tmp = []
        for x in linha.split(','):
            tmp.append(float(x))
        matriz.append(tmp[:])
    arq.close()

    #print(matriz[3][3])
    return matriz

def tamanhoArquivo(nomeArquivo):
    arq = open(nomeArquivo, 'r')
    texto = arq.readlines()
    tamanho=0

    for linha in texto:
        tamanho=tamanho+1
    arq.close()

    return tamanho


#funcao para imprimir uma matriz
def imprimirMatriz(matriz):
    for i in range(len(matriz)):
        for j in range(len(matriz[i])):
            print '  %0.5f' % (matriz[i][j]),

        print '\n'



#resolve a matriz triangular
def resolveTri(matriz):
    linhad = len(matriz) -1
    colunad = len(matriz[0]) -2

    u = linhad
    x = []

    x.append(matriz[linhad][colunad]/(matriz[u][u]))
    u=u-1

    for i in range(u+1):
        linhad=linhad-1
        soma=0
        for j in range(u+1):
            if j<len(x):
                soma=soma+(matriz[u-i][u-j+1]*x[j])
        x.append(((matriz[linhad][colunad])-soma)/matriz[u-i][u-i])
    x.reverse()
    return x

#Funcao que acha o maior numero entre os valores da coluna passada
def achaMaior(mat, i, j):
    maior = 0
    maiorIndice = i
    while i < len(mat):
       # print(mat[i][j])
        if abs(mat[i][j]) > maior:
            maior = mat[i][j]
            maiorIndice = i
        i += 1
    return maiorIndice


#Se pivotacao=0 sem pivotacao pivotacao=1 com pivotacao
def eliminacaoGauss(matriz,pivotacao):
    n=len(matriz) #ordem da matriz

    for i in range(n-1):
        if pivotacao == 1:

            maior=achaMaior(matriz,i,i)
            aux=matriz[i]
            matriz[i]=matriz[maior]
            matriz[maior]=aux

        for j in xrange(i,n-1):
            m=matriz[j+1][i]/matriz[i][i]*(-1)

            for x in xrange(0,n+1):
                aux=(matriz[i][x]*m)+matriz[j+1][x]
                matriz[j+1][x]=aux

    #imprimirMatriz(matriz)
    x=resolveTri(matriz)
    return x

def copiarColuna(matriz, coluna):
    n=len(matriz) #ordem da matriz
    v=[]
    for x in xrange(0,n):
        v.append(matriz[x][coluna])

    return v




def criterioLinhas(matriz):

    u = len(matriz)

    flag=0 #0 - Satisfaz 1 - Nao Satisfaz


    for i in xrange(0,u):
        soma=0
        for j in xrange(0,u):
            if i==j:
                a=matriz[i][j]
            else:
                soma=soma+abs(matriz[i][j])

        #print'soma', (soma), 'a', (a)
        if abs(a)<soma:
            return 1

    return 0
        
        

def jacobi(matriz, xInicial,  toler, iterMax, imprimir):

    aux=criterioLinhas(matriz);
    x=[]
    v=[]

    if aux==1:
        print '\nEsse sistema nao satisfaz o criterio das linhas !!!'
        c = int(input('Deseja continuar mesmo assim? (0-Sim 1-nao)'))

        if c == 1:
            return x, 0

    n=len(matriz) #ordem da matriz
    a=matriz
    b=copiarColuna(matriz, len(matriz[0])-2)
    #print 'b ', (b)
    iter=0
    flag=0

    for i in xrange(0,n):
        r=1/a[i][i]
        for j in xrange(0,n):
            if(i!=j):
                a[i][j]=a[i][j]*r
        b[i]=b[i]*r
        x.append(xInicial)
    while flag == 0:
        iter=iter+1
        for i in xrange(0,n):
            soma=0
            for j in xrange(0,n):
                if i!=j:
                    soma=soma+a[i][j]*x[j]
            v.append(b[i]-soma)
        normaNum=0.0
        normaDen=0.0
        for i in xrange(0,n):
            t=abs(v[i]-x[i])
            if t>normaNum:
                normaNum=t
            if abs(v[i])>normaDen:
                normaDen=abs(v[i])
            x[i]=v[i]
        normaReal=normaNum/normaDen
        if imprimir == 0:
            print '\nIter:', (iter) 
            print '  x:', (x)
            print '  norma:', (normaReal)
        v=[]

        #Testes de convergencia
        if(normaReal<=toler) or (iter>=iterMax):
            flag=1
 
    return (x, iter)


def criterioSassenfeld(matriz):
    u = len(matriz)
    vetorB=[]

    for i in xrange(0,u):
        vetorB.append(1)

    flag=0 #0 - Satisfaz 1 - Nao Satisfaz

    for i in xrange(0,u):
        soma=0
        for j in xrange(0,u):
            if i==j:
                a=abs(matriz[i][j])
            else:
                soma=soma+(abs(matriz[i][j])*vetorB[j])

        b=soma/a
        vetorB[i]=b
       # print'soma', (soma), 'a', (a), 'b', (b)
        if b>1:
            return 1

    return 0



def seidel(matriz,xInicial, toler, iterMax, imprimir):

    x=[]
    v=[]

    aux = criterioSassenfeld(matriz)
    if aux==1:
        print '\nEsse sistema nao satisfaz o criterio de Sassenfeld !!!'
        c = int(input('Deseja continuar mesmo assim? (0-Sim 1-nao)'))

        if c == 1:
            return x, 0 


    n=len(matriz) #ordem da matriz
    a=matriz
    b=copiarColuna(matriz, len(matriz[0])-2)
    #print 'b ', (b)
    iter=0
    flag=0

    for i in xrange(0,n):
        r=1/a[i][i]
        for j in xrange(0,n):
            if(i!=j):
                a[i][j]=a[i][j]*r
        b[i]=b[i]*r
        x.append(xInicial)
    while flag == 0:
        iter=iter+1
        v=[]
        for i in xrange(0,n):
            soma=0
            for j in xrange(0,n):
                if i!=j:
                    soma=soma+a[i][j]*x[j]
            v.append(x[i])
            x[i]=b[i]-soma
        normaNum=0.0
        normaDen=0.0
        for i in xrange(0,n):
            t=abs(x[i]-v[i])
            if t>normaNum:
                normaNum=t
            if abs(x[i])>normaDen:
                normaDen=abs(x[i])
        normaReal=normaNum/normaDen
        if imprimir == 0:
            print '\nIter:', (iter) 
            print '  x:', (x)
            print '  norma:', (normaReal)

        #Testes de convergencia
        if(normaReal<=toler) or (iter>=iterMax):
            flag=1


    return (x,iter)         


def calculaResiduo(matriz, x):
    u = len(matriz)
    b = len(matriz[0])-2
    resp = []


#resolve o sistema
    for i in xrange(0,u):
        soma = 0
        for j in xrange(0,u):
            soma=soma+matriz[i][j]*x[j]
        resp.append(soma)

#Calcula Residuo
    residuo=[]

    for x in xrange(0,u):
        d=matriz[x][b]-resp[x]
        residuo.append(d)


    return residuo

 


def principal():

    ativo=0

    while ativo == 0:

        print'\nMenu:\n'
        print'1-Utilizar algum Metodo'
        print'2-Gerar Graficos'
        print'3-Comparar Residuo'
        print'4-Sair\n'
        menu = int(input('Digite sua opcao: '))

        if menu == 4:
            ativo=1

        if menu == 1:

            print'\nEscolha o Metodo:\n'
            print'1-Eliminacao de Gauss sem pivoteamento'
            print'2-Eliminacao de Gauss com pivoteamento'
            print'3-Gauss-Jacobi'
            print'4-Gauss-Seidel'
            print'5-Voltar\n'

            subMenu = int(input('Digite sua opcao: '))

            if subMenu == 1:

                nomeArquivo = raw_input('Digite o nome do arquivo: ')
                
                matriz = lerArquivo(nomeArquivo)
                print'\nMatriz de Entrada\n'
                imprimirMatriz(matriz)

                print '\n\nEleminacao de Gauss sem Pivotamento'
                resp=eliminacaoGauss(matriz,0)
                imprimirMatriz(matriz)
                print'resposta ',(resp)

            if subMenu == 2:

                nomeArquivo = raw_input('Digite o nome do arquivo: ')
                
                matriz = lerArquivo(nomeArquivo)
                print'\nMatriz de Entrada\n'
                imprimirMatriz(matriz)

                print '\nEleminacao de Gauss com Pivotamento'
                resp=eliminacaoGauss(matriz,1)
                imprimirMatriz(matriz)
                print'resposta ',(resp)

            if subMenu == 3:

                nomeArquivo = raw_input('Digite o nome do arquivo: ')
                erroMax = float(input('Digite o valor do erro maximo: '))
                iterMax = float(raw_input('Digite o valor maximo de iteracoes: '))
                x=0
                
                matriz = lerArquivo(nomeArquivo)
                print'\nMatriz de Entrada\n'
                imprimirMatriz(matriz)


                print '\n\nGauss Jacobi\n'
                r, k=jacobi(matriz,x, erroMax, iterMax,0)
                print '\nResposta Final', (r)
                print 'iteracoes: ', (k)

            if subMenu == 4:

                nomeArquivo = raw_input('Digite o nome do arquivo: ')
                erroMax = float(input('Digite o valor do erro maximo: '))
                iterMax = float(raw_input('Digite o valor maximo de iteracoes: '))
                x=0
                
                matriz = lerArquivo(nomeArquivo)
                print'\nMatriz de Entrada\n'
                imprimirMatriz(matriz)

                print '\n\nGauss Seidel\n'
                r, k =seidel(matriz,x, erroMax, iterMax,0)
                print '\nResposta Final', (r)
                print 'iteracoes: ', (k)


        if menu == 2:

            print'\nEscolha o Metodo:\n'
            print'1-Grafico comparando o tempo'
            print'2-Grafico para metodo iterativo'
            print'3-Voltar\n'

            subMenu = int(input('Digite sua opcao: '))

            if subMenu == 1:

                flag=0
                tamanhoMat=[]
                tempoExecucaoGSP=[]
                tempoExecucaoGCP=[]
                tempoExecucaoGJ=[]
                tempoExecucaoGS=[]


                while flag ==0:

                    nomeArquivo = raw_input('Digite o nome do arquivo: ')

                    matriz = lerArquivo(nomeArquivo)
                    tamanho = tamanhoArquivo(nomeArquivo)
                    tamanhoMat.append(tamanho)

                    tempoInicial = time.time()
                    resp=eliminacaoGauss(matriz,0)
                    tempoFinal = time.time()
                    tempoExecucaoGaussSemPiv=tempoFinal-tempoInicial
                    tempoExecucaoGSP.append(tempoExecucaoGaussSemPiv)


                    matriz = lerArquivo(nomeArquivo)


                    tempoInicial = time.time()
                    resp=eliminacaoGauss(matriz,1)
                    tempoFinal = time.time()
                    tempoExecucaoGausscomPiv=tempoFinal-tempoInicial
                    tempoExecucaoGCP.append(tempoExecucaoGausscomPiv)

                    matriz = lerArquivo(nomeArquivo)

                    x=0
                    erroMax=0.000000000000001
                    iterMax = 50

                    tempoInicial = time.time()
                    r, k =jacobi(matriz,x, erroMax, iterMax,1)
                    tempoFinal = time.time()
                    tempoExecucaofinalGJ=tempoFinal-tempoInicial
                    tempoExecucaoGJ.append(tempoExecucaofinalGJ)



                    matriz = lerArquivo(nomeArquivo)

                    tempoInicial = time.time()
                    r, k =seidel(matriz,x, erroMax, iterMax,1)
                    tempoFinal = time.time()
                    tempoExecucao=tempoFinal-tempoInicial
                    tempoExecucaoGS.append(tempoExecucao)



                    flag = int(input('Deseja adicionar mais algum sistema? (0-Sim , 1-Nao): '))



                plt.plot(tempoExecucaoGCP, tamanhoMat, color='red' , label='Gauss com Pivo')
                plt.plot(tempoExecucaoGSP, tamanhoMat, color='blue' , label='Gauss sem Pivo')
                plt.plot(tempoExecucaoGJ, tamanhoMat, color='green' , label='Gauss Jacobi')
                plt.plot(tempoExecucaoGS, tamanhoMat, color='yellow' , label='Gauss Seidel')
                plt.xlabel("Tempo")
                plt.ylabel("Tamanho")
                plt.title("Grafico de Tempo")
                plt.legend()
                plt.show()

            if subMenu ==2:

                flag=0
                tamanhoMat=[]
                numIteracoesGJ=[]
                numIteracoesGS=[]

                while flag==0:

                    nomeArquivo = raw_input('Digite o nome do arquivo: ')
                    erroMax = float(input('Digite o valor do erro maximo: '))
                    iterMax = float(raw_input('Digite o valor maximo de iteracoes: '))

                    erroVet=[]

                    matriz = lerArquivo(nomeArquivo)
                    tamanho = tamanhoArquivo(nomeArquivo)
                    tamanhoMat.append(tamanho)

                    x=0
                    erroMax=0.000000000000001
                    iterMax = 50
                    r, k=jacobi(matriz, x, erroMax, iterMax,1)
                    numIteracoesGJ.append(k)


                    matriz = lerArquivo(nomeArquivo)

                    x=0
                    r, k =seidel(matriz,x, erroMax, iterMax,1)
                    numIteracoesGS.append(k)


                    flag = int(input('\n\nDeseja adicionar mais algum sistema? (0-Sim , 1-Nao): '))

                plt.plot(numIteracoesGJ, tamanhoMat, color='red' , label='Gauss Jacobi')
                plt.plot(numIteracoesGS, tamanhoMat, color='blue' , label='Gauss Seidel')
                plt.xlabel("Iteracoes")
                plt.ylabel("Tamanho")
                plt.title("Numero de iteracoes x metodo x tamanho do sistema")
                plt.legend()
                plt.show()

        if menu == 3:

            nomeArquivo = raw_input('Digite o nome do arquivo: ')
                
            matriz = lerArquivo(nomeArquivo)
            print'\nMatriz de Entrada\n'
            imprimirMatriz(matriz)

            print '\nEleminacao de Gauss sem Pivotamento'
            resp=eliminacaoGauss(matriz,0)

            matriz = lerArquivo(nomeArquivo)
            residuoGSP=calculaResiduo(matriz, resp)
            print'\nresiduo ',(residuoGSP)

            print '\nEleminacao de Gauss com Pivotamento'
            resp=eliminacaoGauss(matriz,1)

            matriz = lerArquivo(nomeArquivo)
            residuoGCP=calculaResiduo(matriz, resp)
            print'\nresiduo ',(residuoGCP)

            erroMax = 0.000000000000001
            iterMax = 50
            x=0

            print '\nGauss Jacobi'
            r, k=jacobi(matriz,x, erroMax, iterMax,1)

            matriz = lerArquivo(nomeArquivo)
            residuoGJ=calculaResiduo(matriz, r)
            print'\nresiduo ',(residuoGJ)

            print '\nGauss Seidel'
            r, k =seidel(matriz,x, erroMax, iterMax,1)

            matriz = lerArquivo(nomeArquivo)
            residuoGS=calculaResiduo(matriz, r)
            print'\nresiduo ',(residuoGS)

            tamanho = tamanhoArquivo(nomeArquivo)

            rTaxa = [0,0,0,0]

            for x in xrange(0,tamanho):

                if residuoGSP[x]<0:
                    rTaxa[0]=rTaxa[0]+1

                if residuoGCP[x]<0:
                    rTaxa[1]=rTaxa[1]+1

                if residuoGJ[x]<0:
                    rTaxa[2]=rTaxa[2]+1

                if residuoGS[x]<0:
                    rTaxa[3]=rTaxa[3]+1

            maiorTaxa=0
            for x in xrange(0,4):
                if rTaxa[x]>maiorTaxa:
                    iMaior=x
                    maiorTaxa=rTaxa[x]

            if iMaior == -1:
                print '\n\nOs Metodos tem os mesmos residuos!!!!\n'
            else:
                if iMaior==0:
                    print '\n\nGauss sem pivotamento tem mais residuos!!!!\n'

                if iMaior==1:
                    print '\n\nGauss com pivotamento tem mais residuos!!!!\n'

                if iMaior==2:
                    print '\n\nGauss Jacobi tem mais residuos!!!!\n'

                if iMaior==3:
                    print '\n\nGauss Seidel tem mais residuos!!!!\n'
        





principal()
