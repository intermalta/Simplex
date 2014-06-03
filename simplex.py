#!/usr/bin/python

from numpy import *
import scipy as Sci
import scipy.linalg

def pivot(i,j, base, nbase):
	aux = base[i]
	#print aux
	base[i] = nbase[j]
	print "variavel %d entra na base" %j
	#print base
	nbase[j] = aux
	print "variavel %d sai da base" %aux
	#print nbase
	return base, nbase

def simplex(A, b, c):
	[m,n] = A.shape
	#print [m,n]

	Aext = concatenate((A, eye(m)), 1)
	#print Aext
	
	cext = concatenate((c.transpose(), zeros(m)))
	#print cext
	#print c

	bmin = min(b)
	#print bmin
	if bmin < 0:
		#implementar fase 1 do simplex
		print "Implementar fase 1"
	else:
		base = array([i for i in range(n, n+m)])
		#print base
		nbase = array([i for i in range(0, n)])
		#print nbase
	
	#simplex fase 2
	flag = 0
	while flag == 0:
		
		B = Aext[:,base]
		#print B
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
			flag = 1
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
			else:
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
				posicaoLinha = nanargmin(temp)
				if posicaoLinha == -1:
					print "Menor elemento nao encontrado"
					exit(-2)
				base, nbase = pivot(posicaoLinha, posicaoColuna, base, nbase)
if __name__ == "__main__":
	#c = array([4, 3])
	#A = array([[2, 1], [1,2]])
	#b = array([4,4])

	c = array([4,3])
	A = array([[1,-1], [2,-1], [0,1]])
	b = array([1,3,5])
		
	xb, z, y =simplex(A, b, c)
	print xb
	print z
	print y


