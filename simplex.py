#!/usr/bin/python

from numpy import *
from random import randint
import time


eps = 0.000001 #tolerancia utilizada para comparacao de precisao

def print_tableau(B, NonB, bi, base, nbase):
    Binv = B.copy()
    N = NonB.copy()
    b = bi.copy()
    tableau = zeros(shape=(B.shape[0], B.shape[1]+N.shape[1]), dtype=float)
    tableau[:,base] = eye(B.shape[0])
    tableau[:,nbase] = dot(Binv,N)
    print tableau

def getTableau(Binv, N, b, base, nbase):
    tableau = zeros(shape=(Binv.shape[0], Binv.shape[1]+N.shape[1]), dtype=float)
    tableau[:,base] = eye(Binv.shape[0])
    tableau[:,nbase] = dot(Binv,N)
    return tableau

def randomGraph():
    nodes = randint(5,25)
    print 'Nodes: ' +str(nodes)
    edges = 0
    graph = zeros(shape=(nodes,nodes), dtype=float)
    for i in xrange(nodes):
        for j in xrange(nodes):
            if i != j: #do not creat loops
                e = randint(0,1)
                graph[i][j] = e
                if (e == 1):
                    edges = edges + 1
                
    print 'Edges total= %d' %edges
                
    return graph, nodes, edges

def getRestrictionsPLI():
    g,n,m = randomGraph()
    A = zeros(shape=(m,n), dtype=float)
    row = 0
    for i in xrange(n):
        for j in xrange(n):
            if g[i][j] == 1:
                A[row][i] = 1
                A[row][j] = 1
                row = row + 1
    b = ones(m)*(-1)
    c = ones(n)*(-1) # adaptando um problema de minimizacao a um problema de maximizacao
    return A*(-1),b,c

def validate_input(A, b, c):
    [m, n] = A.shape
    if (m != len(b)):
        return 1
    if (n != len(c)):
        return 1
    return 0

def pivot(i, j, base, nbase):
   
    #print 'base antes do pivoteamento -> ' + (base+1).__str__()
    
    b = base.copy()
    n = nbase.copy()
    
    aux = b[i]
    b[i] = n[j]
    n[j] = aux
    
    #print 'base depois do pivoteamento -> ' + (b+1).__str__()
    return b, n

def revised_row(B, u, l):
    [m, n] = B.shape
    if l > n:
        print 'error with row limit in revised_row'
        exit()
    if (u.shape[0] < l):
        print 'error with vector dimension'
    Binv = B.copy()
    for i in range(m):
        for j in range(n):
            if i != l:
                Binv[i][j] = B[i][j] - ((B[l][j] * u[i]) / u[l])
            else:
                Binv[i][j] = B[l][j] / u[l] 

    return Binv
    
def get_var_sainte(direcao, Binv, b):
    u = direcao.copy()
        
    for i in range(len(u)):
        if u[i] <= 0:
            u[i] = nan

    temp = dot(Binv, b)
    temp = temp / u
    
    return nanargmin(temp)

def phase_one(Aext, b, bmin, m, n):
        # adding w to A matrix at last column
        w = ones((m, 1)) * (bmin)
        
        AwithW = concatenate((Aext, w), 1)
        # print AwithW
        
        #add -1 to w position in the new cost vector
        cwithW = concatenate((zeros(m + n), ones(1) * (-1)))
        base1 = array([i for i in range(n, n + m)])
        nbase1 = array([i for i in range(0, n + 1)])
        #colocando a variavel w no vetor de variaveis nao basicas
        nbase1[-1] = n + m
        posW = n + m
        entrante = len(nbase1) -1

        # adding w to the base
        sainte = argmin(b)
        base1, nbase1 = pivot(sainte, entrante, base1, nbase1)

        Bwinv = revised_row(eye(m), w, sainte)
        tentativas = 200
        while tentativas > 0:
            tentativas = tentativas - 1
            Nw = AwithW[:, nbase1]

            cbw = cwithW[:, base1]
            cnw = cwithW[:, nbase1]

            yw = dot(cbw, Bwinv)
            credw = cnw - dot(yw, Nw)
            xbw = dot(Bwinv, b)
            zw = dot(cbw, xbw)
            
            if max(credw) < 0 or zw > 0:
                print 'Problema inviavel'
                return(-1, xbw)
            
            if zw == 0:
                # we are at optimal z -> search for column w
                posw = -1
                for i in (range(len(nbase1))):
                    if nbase1[i] == posW:
                        posw = i
                        break
                if posw == -1:
                    #base degenerada
                    #procurando w na base
                    for i in (xrange(len(base1))):
                        if base1[i] == posW:
                            posw = i
                            break
                    #atualizando B-1
                    base1, nbase1 = pivot(posw, posw, base1, nbase1)
                    Bwinv = revised_row(Bwinv, dot(Bwinv, Nw[:,posw]), posw) 
                    print 'w not found in non bases vector'
                    #exit(3)

                Nw = delete(Nw, posw)
                # print Nw
                nbase1 = delete(nbase1, posw)
                print '====='
                print 'initial base'
                print base1
                print 'Fim fase 1 com %d iteracoes' % (200 - tentativas)
                return 3,base1, nbase1, Bwinv
                    
            posicaoColunaw = -1
            for j in range(len(credw)):
                if credw[j] > 0:
                    posicaoColunaw = j
                    break
            
            # direcao do deslocamento
            direcaow = dot(Bwinv, Nw[:, posicaoColunaw])
            dwmax = direcaow.max()
        
            if dwmax < 0:
                print "Problema Ilimitado"
            else:
                posicaoLinhaw = get_var_sainte(direcaow, Bwinv, b)
                base1, nbase1 = pivot(posicaoLinhaw, posicaoColunaw, base1, nbase1)
            Bwinv = revised_row(Bwinv, direcaow, posicaoLinhaw)
        
        #maximo de iteracoes atingido
        return(-2, xbw)
             
def simplex(A, b, c):
    
    if validate_input(A, b, c):
        print 'input with bounds problems'
        exit(-1)
    [m, n] = A.shape
    # print [m,n]
    Aext = concatenate((A, eye(m)), 1)
    # print Aext
    cext = concatenate((c.transpose(), zeros(m)))

    bmin = min(b)
    # print bmin
    if bmin < 0:
        print "phase 1"
        ret = phase_one(Aext, b, bmin, m, n)
        if ret[0] == -1:
            return ret
        elif ret[0] == -2 or ret[0] == 0:
            print ret 
            exit(-3)
        else:
            base = ret[1]
            nbase = ret[2]
            Binv = ret[3]

    else:
        base = array([i for i in range(n, n + m)])
        nbase = array([i for i in range(0, n)])
        B = Aext[:, base]
        Binv = B.copy()
        
        
    # print Binv

    # simplex fase 2
    print 'inicio fase 2'
    tentativas = 200
    while tentativas > 0:
        tentativas = tentativas -1
        N = Aext[:, nbase]
        cb = cext[:, base]
        cn = cext[:, nbase]
        y = dot(cb, Binv)
        # print y
        cred = cn - dot(y, N)
        # print cred
        vec = cred[cred > 0]
        xb = dot(Binv, b)
        z = dot(cb, xb)

        if vec.sum() == 0:
            print base, nbase
            print 'Fim fase 2 - %d iteracoes' % (200 - tentativas)
            #t = getTableau(Binv,N,b, base, nbase)
            return 1,xb, z, y, cred, base, nbase#, t, dot(Binv,b)

        # anti cycle - Bland's rule
        posicaoColuna = -1
        #posicaoColuna = argmax(cred)
        for j in range(len(cred)):
            if cred[j] > 0:
                maximo = cred[j]
                posicaoColuna = j
                break
        
        # direcao do deslocamento
        direcao = dot(Binv, N[:, posicaoColuna])
        dmax = direcao.max()
    
        if dmax < 0:
            print "Problema Ilimitado"
            return (0, direcao*(-1))
            
        else:
            posicaoLinha = get_var_sainte(direcao, Binv, b)
            base, nbase = pivot(posicaoLinha, posicaoColuna, base, nbase)
        #print 'base: ' + (base+1).__str__()
        #print 'nbase' + (nbase+1).__str__()
        #print_tableau(Binv, N, b, base, nbase)
        Binv = revised_row(Binv, direcao, posicaoLinha)
           
    #maximo de iteracoes atingido
    return (-2, xb)

def print_relatorio(ret,A):
    if (ret[0] == -2):
        print 'Estado: '+ str(ret[0])
        print'Maximo de iteracoes atingindo'
        print 'Ultima solucao obtida: ' + ret[1].__str__()
    elif (ret[0] == -1):
        print 'Estado: '+ str(ret[0])
        print 'Problema Inviavel'
        print 'Ultima solucao obtida: ' + ret[1].__str__()
    elif (ret[0] == 0):
        print 'Estado: '+ str(ret[0])
        print 'Problema ilimitado'
        print 'Direcao Encontrada: ' + ret[1].__str__()
    elif (ret[0] == 1):
        print 'Estado '+ str(ret[0])
        print 'valor otimo: ' + ret[2].__str__()
        #print 'xb: ' + ret[1].__str__()
        base = array([0. for i in range(A.shape[1])]) #n
        nbase = array([0. for i in range(A.shape[0])]) #m
        for i in range(len(ret[5])):
            if ret[5][i] > A.shape[1] -1:
                nbase[ret[5][i] - A.shape[0]] = ret[1][i]
            else:
                base[ret[5][i]] = ret[1][i]
        print 'variaveis: ' + (base).__str__()
        print 'variaveis de folga' + (nbase).__str__()    
        print 'custos reduzidos ' + ret[4].__str__()
        print 'variaveis duais: ' + ret[3].__str__()
        print 'indices de variaveis basicas: ' + (ret[5]+1).__str__()
        print 'indices de variaveis nao basicas: ' + (ret[6]+1).__str__()
        
def get_xb_var_originais(xb, base, n):
    
    xbOriginal = array([0. for i in xrange(n)])
    
    for i in xrange(len(base)):
        if base[i] < n:
            xbOriginal[base[i]] = xb[i]

    return xbOriginal

#testa o valor de B^(-1)*b para a verificacao da restricao dos numeros serem inteiros   
def testa_inteiros(b):
    for i in b:
        if abs(i - round(i)) > eps:
            return False
    return True
    
def getLinhaGomory(t, indice):
    l = array([0. for i in xrange(len(t))]) 
    
    #for i in xrange(len(t)):
    if t[indice] > 0:
        l[indice] = 1
        value = math.floor(t[indice])
            #break
    
    #value = rhs - math.floor(rhs)
    
    return l, value
    
def getCut(t):
    
    rhs = array([])
    for i in xrange(len(t)):
        if (abs(t[i] - round(t[i])) >eps):
            cut,value = getLinhaGomory(t, i)
            rhs = append(rhs, float(value)) #*(-1))
            break     
    

    #multiplicando o corte por -1 para transformar a restricao do tipo <=
    return cut,rhs
            
            
            
if __name__ == "__main__":
    #c = array([4, 3])
    #A = array([[2, 1], [1,2]])
    #b = array([4,4])

    # c = array([4,3])
    # A = array([[1,-1], [2,-1], [0,1]])
    # b = array([1,3,5])
        
    #c = array([4, 3])
    #A = array([[2, 1], [1, 2], [-1, -1]])  
    #b = array([4, 4, -1])
    
    
    #########Chvatal pg 39
    #c = array([1,-1,1])
    #A = array([[2,-1,2],[2,-3,1], [-1,1,-2]])
    #b = array([4,-5,-1])
    #ret = simplex(A,b,c)
    #print_relatorio(ret,A)
  
    
    ############Exemplo de problema ilimitado
    #c = array([2.,3.])
    #A = array([[-1.,2.],[-5.,1.]])
    #b = array([4.,1.])
    #ret = simplex(A,b,c)
    #print_relatorio(ret)
    
    ##### GOMORY
      
    #c = array([4,-1])
    #A = array([[7,-2],[0,1], [2,-2]])
    #b = array([14,3,3])
    #while(True):
    #testes grafo
    A,b,c = getRestrictionsPLI()
    start = time.clock()
    ret = simplex(A,b,c)
    print_relatorio(ret,A)
    print 'time elapsed: ' + str(time.clock() - start) 
    #A = concatenate((A, eye(A.shape[0])),1)
    estado = ret[0]
    maximo = 0
    if estado == 1:
        
        #testa se os valores das variaveis sao inteiras
        while(testa_inteiros(ret[1]) == False and maximo < 10):
            start_gomory = time.clock()
            t = get_xb_var_originais(ret[1], ret[5], A.shape[1])
            cut, rhs = getCut(t)
        
            A = vstack((A, cut))
            b = concatenate((b, rhs))
            ret = simplex(A,b,c)
            maximo = maximo + 1
        
        print_relatorio(ret,A)
        print 'Elapsed gomory: ' + str(time.clock() - start_gomory)
        del(A)
        del(b)
        del(c)
    # #teste matriz revisada
    # B = array([[1., 2., 3.], [-2., 3., 1.], [4.,-3.,-2.]])
    # u = array([-4,2,2])
    # revised_row(B,u,2)
    # B = array([[1., 0., 0.], [2., 1., 0.], [0., 0., 1.]])
    # u = array([1,2,0])

    # print revised_row(B, u, 0)
    # print linalg.inv(B)


    print '==== resposta ===='
    #print xb
    #print z
    #print y
    
    
    