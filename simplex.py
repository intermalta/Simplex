#!/usr/bin/python

from numpy import *
from random import randint

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
    nodes = randint(5,10)
    print nodes
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
    for i in xrange(m):
        for j in xrange(m):
            if g[i][j] == 1:
                A[row][i] = 1
                A[row][j] = 1
                row = row + 1

def validate_input(A, b, c):
    [m, n] = A.shape
    if (m != len(b)):
        return 1
    if (n != len(c)):
        return 1
    return 0

def pivot(i, j, base, nbase):
   
    print 'base antes do pivoteamento -> ' + (base+1).__str__()
    b = base.copy()
    n = nbase.copy()
    aux = b[i]
    b[i] = n[j]
    n[j] = aux
    print 'base depois do pivoteamento -> ' + (b+1).__str__()
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

        # add -1 to w position in the new cost vector
        cwithW = concatenate((zeros(m + n), ones(1) * (-1)))

        base1 = array([i for i in range(n, n + m)])
        
        nbase1 = array([i for i in range(0, n + 1)])
        #colocando a variavel w no vetor de variaveis nao basicas
        nbase1[-1] = n + m

        # adding w to the base
        entrante = len(nbase1) -1 
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
            
            if max(cbw) < 0 or zw > 0:
                print 'Problema inviavel'
                return(-1, xbw)
            
            if zw == 0:
                # we are at optimal z -> search for column w
                posw = -1
                for i in (range(len(nbase1))):
                    if nbase1[i] == n + m:
                        posw = i
                if posw == -1:
                    print 'w not found in non bases vector'
                    exit(3)

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
            if posicaoColunaw == -1:
                print 'nao achou uma coluna para entrar'
                exit(-20)
            # direcao do deslocamento
            direcaow = dot(Bwinv, Nw[:, posicaoColunaw])
            dwmax = direcaow.max()
        
            if dwmax < 0:
                print "Problema Ilimitado"
            else:
                posicaoLinhaw = get_var_sainte(direcaow, Bwinv, b)
                base1, nbase1 = pivot(posicaoLinhaw, posicaoColunaw, base1, nbase1)
            Bwinv = revised_row(Bwinv, direcaow, posicaoLinhaw)
             
def simplex(A, b, c, gomory=0, cortes=0):
    
    if validate_input(A, b, c):
        print 'input with bounds problems'
        exit(-1)
    [m, n] = A.shape
    # print [m,n]
    if gomory == 0:
        Aext = concatenate((A, eye(m)), 1)
        # print Aext
        cext = concatenate((c.transpose(), zeros(m)))
    else:
        if cortes == 0:
            print 'Iteracao sem corte novo'
            exit()
        Aext = concatenate((A, eye(m,cortes)*(1)),1)
        cext = concatenate((c, zeros(cortes)))

    
    bmin = min(b)
    # print bmin
    if bmin < 0:
        print "phase 1"
        ret = phase_one(Aext, b, bmin, m, n)
        if ret[0] == -1:
            return ret
        else:
            base = ret[1]
            nbase = ret[2]
            Binv = ret[3]

    else:
        if gomory == 0:
            base = array([i for i in range(n, n + m)])
            nbase = array([i for i in range(0, n)])
        else:
            base = array([i for i in range(n, n + cortes)])
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
            t = getTableau(Binv,N,b, base, nbase)
            return 1,xb, z, y, cred, base, nbase, t, dot(Binv,b)

        # anti cycle - Bland's rule
        posicaoColuna = -1
        for j in range(len(cred)):
            if cred[j] > 0:
                #maximo = cred[j]
                posicaoColuna = j
                break
        if posicaoColuna == -1:
            print 'Erro ao achar a variavel que entra na base'
            exit(-3)
        # direcao do deslocamento
        direcao = dot(Binv, N[:, posicaoColuna])
        dmax = direcao.max()
    
        if dmax < 0:
            print "Problema Ilimitado"
            return (0, direcao*(-1))
            
        else:
            posicaoLinha = get_var_sainte(direcao, Binv, b)
            base, nbase = pivot(posicaoLinha, posicaoColuna, base, nbase)
        
        Binv = revised_row(Binv, direcao, posicaoLinha)   
    #maximo de iteracoes atingido
    return (-2, xb)

def print_relatorio(ret):
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
            if ret[5][i] > A.shape[0]:
                nbase[ret[5][i] - A.shape[0]] = ret[1][i]
            else:
                base[ret[5][i]] = ret[1][i]
        print 'variaveis: ' + (base).__str__()
        print 'variaveis de folga' + (nbase).__str__()    
        print 'custos reduzidos ' + ret[4].__str__()
        print 'variaveis duais: ' + ret[3].__str__()
        print 'indices de variaveis basicas: ' + (ret[5]+1).__str__()
        print 'indices de variaveis nao basicas: ' + (ret[6]+1).__str__()
   
def testa_inteiros(b, eps):
    for i in b:
        if abs(i - round(i)) > eps:
            return False
    return True
    
def getLinhaGomory(t, rhs):
    l = t.copy() 
    
    for i in xrange(len(t)):
        l[i] = t[i] - math.floor(t[i])
    
    value = rhs - math.floor(rhs)
    
    return l, value
    
def getCuts(t, b, eps):
    cut = array([])
    rhs = array([])
    primeiro = True
    slack = 0
    
    for i in xrange(len(b)):
        #if (abs(b[i] - round(b[i])) > eps):
        l, value = getLinhaGomory(t[i], b[i])
        if (primeiro == True): 
            cut = concatenate((cut, l))
            primeiro = False
        else:
            cut = concatenate((cut, l)).reshape(slack + 1, t.shape[1])

        rhs = append(rhs, float(value))
        slack = slack + 1
            
    #cut = concatenate((cut, eye(slack)), 1)
    return cut,rhs
            
            
            
if __name__ == "__main__":
    # c = array([4, 3])
    # A = array([[2, 1], [1,2]])
    # b = array([4,4])

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
  
    
    ############Exemplo de problema ilimitado
    #c = array([4,3])
    #A = array([[-2,-1],[-1,-2]])
    #b = array([-4,-4])
    #c = array([1,3])
    #A = array([[-1,-1],[-1,1],[-1,2]])
    #b = array([-3,-1,2])
    #c = array([10,20])
    #A = array([[2,3],[4,1]])
    #b = array([6,4])
    #ret = simplex(A,b,c)
    #print_relatorio(ret)
    
    
    
    #A,m = randomGraph()
    #print A 
    #xb, z, y = simplex(A, b, c)
    ##### GOMORY  
    c = array([5,8])
    A = array([[1,1],[5,9]])
    b = array([6,45])
    #while(True):
        
    eps = 0.00001
    ret = simplex(A,b,c)
    print_relatorio(ret)
    A = concatenate((A, eye(A.shape[0])),1)
    
    #testa se os valores das variaveis sao inteiras
    while(testa_inteiros(ret[8], eps) == False):
        cut, rhs = getCuts(ret[7], ret[8], eps)
        [m,n] = cut.shape
        #A = concatenate((A, zeros((A.shape[0], m))), 1)
        A = concatenate((cut, A))
        
        b = concatenate((b, rhs))
        c = concatenate((c, zeros(A.shape[1] - len(c))))
        ret = simplex(A,b,c, 1, m)
    
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
