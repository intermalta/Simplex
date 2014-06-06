#!/usr/bin/python

from numpy import *
import scipy as Sci
import scipy.linalg

def print_dictionary(Binv, N, b, cb, cn, xn):
	return
def validate_input(A,b,c):
	[m,n] = A.shape
	if (m != len(b)):
		return 1
	if (n != len(c)):
		return 1
	return 0

def pivot(i,j, base, nbase):
	aux = base[i]
	#print aux
	#print "variavel %d sai da base" %(base[i])
	base[i] = nbase[j]
	#print "variavel %d entra na base" %(nbase[j])
	#print base
	nbase[j] = aux
	#print nbase
	return base, nbase

def revised_row(B, u, l):
	[m,n] =B.shape
	if l > n:
	  print 'error with row limit in revised_row'
	  exit()
	if (u.shape[0] < l):
		print 'error with vector dimension'
	Binv = B
	for i in range(m):
		for j in range(n):
			if i != l:
			  Binv[i][j] = B[i][j] - ((B[l][j]*u[i])/u[l])
			else:
			  Binv[i][j] = B[l][j]/u[l] 

	return Binv
	
def get_var_sainte(direcao, Binv, b):

		for i in range(len(direcao)):
			if direcao[i] <= 0:
				direcao[i] = nan
		#print direcao
		#print Binv
		#print b
		temp = dot(Binv, b)
		#print temp
		temp = temp/direcao
		#print temp
		return nanargmin(temp)

def phase_one(Aext, b, bmin, cext, m,n):
		#adding w to A matrix at last column 
		w = ones((m,1))*(bmin)
		AwithW = concatenate((Aext, w), 1)
		#print AwithW

		#add -1 to w position in the new cost vector
		cwithW = concatenate((zeros(m+n), ones(1)*bmin))
		#print cwithW

		base1 = array([i for i in range(n, n+m)])
		#print base1
		nbase1 = array([i for i in range(0,n+1)])
		nbase1[-1] = n+m
		#print nbase1

		#adding w to the base
		entrante = len(base1) -1  
		#print 'entrante %d' %entrante
		sainte = argmin(b)
		#print 'sainte %d' %sainte
		base1, nbase1 = pivot(sainte, entrante, base1, nbase1)
		print base1, nbase1

		tentativas = 200
		while tentativas > 0:
			Bw = AwithW[:,base1]
			print Bw
			Nw = AwithW[:,nbase1]
			print Nw

			cbw = cwithW[:,base1]
			cnw = cwithW[:,nbase1]

			Bwinv =  linalg.inv(Bw)
			yw = dot(cbw, Bwinv)
			credw = cnw - dot(yw, b)
			xbw = dot(Bwinv, b)
			zw = dot(cbw, xbw)
			
			if max(cbw) < 0 or zw > 0:
				print 'Problema inviavel'
				exit()
			else:
				if zw == 0:
				  #we are at optimal z -> search for column w
				  posw = -1
				  for i in (range(len(nbase1))):
					  if nbase1[i] == n+m:
						  posw = i
				  if posw == -1:
					  print 'w not found in non bases vector'
					  exit(3)

				  Nw = delete(Nw, posw)
				  #print Nw
				  nbase1 = delete(nbase1, posw)
				  print '====='
				  print 'initial base'
				  print base1
				  print 'Fim fase 1 com %d iteracoes' %(200 - tentativas +1)
				  return base1, nbase1
			
			maximow = 0
			for j in range(credw.shape[0]):
				if credw[j] > maximow:
					maximow = credw[j]
					posicaoColunaw = j
					break
		
			#direcao do deslocamento
			direcaow = dot(Bwinv, Nw[:,posicaoColunaw])

			#print "=========="
			#print direcaow
			dwmax = direcaow.max()
		
			if dwmax < 0:
				print "Problema Ilimitado"
			else:
				posicaoLinhaw = get_var_sainte(direcaow, Bwinv, b)
				base1, nbase1 = pivot(posicaoLinhaw, posicaoColunaw, base1, nbase1)
			tentativas = tentativas - 1
def simplex(A, b, c):
	
	if validate_input(A,b,c):
		print 'input with bounds problemns'
		exit(-1)
	[m,n] = A.shape
	#print [m,n]

	Aext = concatenate((A, eye(m)), 1)
	#print Aext
	
	cext = concatenate((c.transpose(), zeros(m)))
	#print cext
	#print c
	#print '================='

	bmin = min(b)
	#print bmin
	if bmin < 0:
		print "phase 1"
		base, nbase = phase_one(Aext, b, bmin, cext, m, n)

	else:
		base = array([i for i in range(n, n+m)])
		#print base
		nbase = array([i for i in range(0, n)])
		#print nbase
	
	
	#simplex fase 2
	print 'inicio fase 2'
	tentativas = 200
	while tentativas > 0:
		B = Aext[:,base]
		print B
		N = Aext[:,nbase]
		#print N
		
		cb = cext[:,base]
		#print cb
		cn = cext[:,nbase]
		#print cn
		
		###### MUDAR AS INVERSOES DE MATRIZES #######
		Binv = linalg.inv(B)
		#print Binv

		y = dot(cb,Binv)
		#print y
		cred = cn - dot(y,N)
		#print cred
		vec=cred[cred > 0]

		if vec.sum() == 0:
			#estamos no otimo
			xb = dot(Binv,b)
			z = dot(cb,xb)
			print base, nbase
			print 'Fim fase 2 - %d iteracoes' %(200 - tentativas +1)
			return xb,z,y
		else:

			#anti cycle - lexicographic order
			maximo = 0
			for j in range(cred.shape[0]):
				if cred[j] > maximo:
					maximo = cred[j]
					posicaoColuna = j
					break
		
			#direcao do deslocamento
			direcao = dot(Binv, N[:,posicaoColuna])
	
			#print "=========="
			#print direcao
			dmax = direcao.max()
		
			if dmax < 0:
				print "Problema Ilimitado"
				exit()
			else:
				posicaoLinha = get_var_sainte(direcao, Binv, b)
				base, nbase = pivot(posicaoLinha, posicaoColuna, base, nbase)
			
			tentativas = tentativas - 1
if __name__ == "__main__":
	 #c = array([4, 3])
	 #A = array([[2, 1], [1,2]])
	 #b = array([4,4])

	 #c = array([4,3])
	 #A = array([[1,-1], [2,-1], [0,1]])
	 #b = array([1,3,5])
		
	c = array([4, 3])
	A = array([[2, 1], [1,2], [-1,-1]])  
	b = array([4,4, -1])
	xb, z, y =simplex(A, b, c)
  
	##teste matriz revisada
	#B = array([[1., 2., 3.], [-2., 3., 1.], [4.,-3.,-2.]])
	#u = array([-4,2,2])
	#revised_row(B,u,2)
	#B = array([[1., 0., 0.], [2., 1., 0.], [0., 0., 1.]])
	#u = array([1,2,0])

	 #print revised_row(B, u, 0)
	 #print linalg.inv(B)


	print '==== resposta ===='
	print xb
	print z
	print y


