#!/usr/bin/python

from numpy import *
import scipy as Sci
import scipy.linalg

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
                #### XXX - excluir a coluna relativa a w
                exit()
            
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
                return base1, nbase1, Bwinv
                    
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
            
def simplex(A, b, c):
    
    if validate_input(A, b, c):
        print 'input with bounds problemns'
        exit(-1)
    [m, n] = A.shape
    # print [m,n]

    Aext = concatenate((A, eye(m)), 1)
    # print Aext
    
    cext = concatenate((c.transpose(), zeros(m)))
    # print cext
    # print c
    # print '================='

    bmin = min(b)
    # print bmin
    if bmin < 0:
        print "phase 1"
        first_fase = True
        base, nbase, Binv = phase_one(Aext, b, bmin, m, n)

    else:
        base = array([i for i in range(n, n + m)])
        # print base
        nbase = array([i for i in range(0, n)])
        # print nbase
    if first_fase == False:
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

        if vec.sum() == 0:
            # estamos no otimo
            xb = dot(Binv, b)
            z = dot(cb, xb)
            print base, nbase
            print 'Fim fase 2 - %d iteracoes' % (200 - tentativas)
            return xb, z, y

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
            exit()
        else:
            posicaoLinha = get_var_sainte(direcao, Binv, b)
            base, nbase = pivot(posicaoLinha, posicaoColuna, base, nbase)
        
        Binv = revised_row(Binv, direcao, posicaoLinha)
            
          
if __name__ == "__main__":
    # c = array([4, 3])
    # A = array([[2, 1], [1,2]])
    # b = array([4,4])

    # c = array([4,3])
    # A = array([[1,-1], [2,-1], [0,1]])
    # b = array([1,3,5])
        
    c = array([4, 3])
    A = array([[2, 1], [1, 2], [-1, -1]])  
    b = array([4, 4, -1])
    
    
    #########Chvatal pg 39
    #c = array([1,-1,1])
    #A = array([[2,-1,2],[2,-3,1], [-1,1,-2]])
    #b = array([4,-5,-1])
    
    xb, z, y = simplex(A, b, c)
  
    # #teste matriz revisada
    # B = array([[1., 2., 3.], [-2., 3., 1.], [4.,-3.,-2.]])
    # u = array([-4,2,2])
    # revised_row(B,u,2)
    # B = array([[1., 0., 0.], [2., 1., 0.], [0., 0., 1.]])
    # u = array([1,2,0])

    # print revised_row(B, u, 0)
    # print linalg.inv(B)


    print '==== resposta ===='
    print xb
    print z
    print y
