from sympy import *
from copy import copy, deepcopy
# para python 3




def printmatrix(a):
    # a es matriz de sympy
    # imprime a por filas
    data = []
    n = a.shape[0] # filas
    m = a.shape[1] # columnas
    for i in range(n):
        data.append([])
        for j in range(m):
            data[i].append(a[m * i + j])
            # print a[m*i + j],' ',
            # print ''
    col_width = max(len(str(word)) for row in data for word in row) + 2  # padding
    for row in data:
        print ('[  ',"".join(str(word).ljust(col_width) for word in row),']')

def matrix(a):
	# a es matriz de sympy
	#formatea por filas
	data = []
	n = a.shape[0] # filas
	m = a.shape[1] # columnas
	for i in range(n):
		data.append([])
		for j in range(m):
			data[i].append(a[m * i + j])
			# print a[m*i + j],' ',
			# print ''
	col_width = max(len(str(word)) for row in data for word in row) + 2  # padding
	ret = ''
	for row in data:
		ret = ret +  '[  '
		ret = ret + "".join(str(word).ljust(col_width) for word in row)
		ret = ret + ']\n'
	return ret[:-1]


def dirac(i,j):
	ret = 0
	if i ==j:
		ret =1
	return ret
	
def column2vector(u):
	# pre: ingresa una matriz de una columna Matrix([[a1],[a2],...,[an]]) 
	# post: devuelve [a1,a2,...,an]
	x = u.T.tolist()
	return x
def columns2matrix(w):
	# pre: w una lista matrices columna,
	# post: devuleve la matriz formada por esas columnas
	c = []
	for i in range(len(w)):
		c.append(w[i].T.tolist()[0])
	r = Matrix(c).T
	return r
	
# ====== INICIO: Operaciones elementales de matrices ================

def E(i,r,j,A):
	# fila_i(A) = fila_i(A) + r*fila_j(A)
	R = []
	rFj = list(r*A.row(j))
	for k in range(i):
		R.append(list(A.row(k)))
	R.append(list(A.row(i) + r*A.row(j)))
	for k in range(i+1,A.shape[0]):
		R.append(list(A.row(k)))
	return Matrix(R)
	
def F(r,j,A):
	# fila_j(A) = r*fila_j(A)
	return E(j,r-1,j,A)
	
def G(i,j,A):
	# fila_i(A) = fila_j(A), fila_j(A) = fila_i(A)
	R = E(i,1,j, A)
	R = E(j,-1,i,R)
	R = E(j,-2,j,R)
	R= E(i,-1,j,R)
	return R

# ====== FIN: Operaciones elementales de matrices ================


def teoremaVS(A):
	# pre: A  es una matriz real
	# post: devuelve [[v0,..,v(n-1)],[w(0),...,w(m-1)], [a0,...,a(r-1)]] tal que
	#	1) v es BON 
	#	2) w es BON
	#	3) ai >0 (0 <= i < r)
	#	4) A(vi) = ai*wi (0 <= i < r), A(vi) = 0 (i >= r)
	n = A.shape[0] # filas
	m = A.shape[1] # columnas
	ATA = A.T*A
	a = [] # lista de autovalores de A^*A
	v0 = [] # lista de autovectores de A^*A
	Ev = ATA.eigenvects() #Ev[0] = autovalor, Ev[1] = multiplicidad, Ev[2] = lista de autov3ectores
	# los vectores se representan como matrices columna. 
	Ev.sort()
	Ev.reverse()
	
	rg = min(n,m) # rango de ATA
	for i in range(len(Ev)):
		for j in range(Ev[i][1]):
			a.append(Ev[i][0])
			v0.append(Ev[i][2][j])
			if Ev[i][0] == 0:
				rg = min(len(v0)-1,rg)
	v0 = GramSchmidt(v0, True)
	
	# Autovalores y autovectores de A^*A: "multiplicidad autovector"'
	v = []
	for i  in range(len(a)):
		v.append(v0[i].T.tolist()[0])
	# print  a
	# print v
	# print rg

	w0 = []
	for i in range(rg):
		w0.append(A*v0[i])
	# print 'w000000000000',w0
	w0 = GramSchmidt(w0, True)
	#for i in range(len(w0)):
	#	print 'w['+str(i)+'] =',w0[i].T.tolist()[0]
	
	# '\nCompletar base w a BON'
	e = []
	for  k in range(n):
		e.append(Matrix(n, 1, lambda i,j: dirac(i,k)) )# base canonica
		w0.append(e[k])
	A1 =  columns2matrix(w0)
	#print 'A1.T', A1.T

	mat = A1
	li = mat.rref()[1] # indice de los vectores que son li
	w1 = []
	for i in li:
		w1.append(w0[i])
	# w1 es una base que completa w
		
	w2 = GramSchmidt(w1, True) # base completada y luego Gram-Scmidt
	w = []
	for i in range(len(w2)):
		w.append(w2[i].T.tolist()[0])
	
	return [v,w,a]
e = [] 
e.append(Matrix([[1],  [0],[0]]))
e.append(Matrix([[0],  [1],[0]]))
e.append(Matrix([[0],  [0],[1]]))

a1, a2, a3 = symbols('a1, a2, a3')
t =  symbols('t')



"""

print '\n\nDescomposicion polar'
A  = Matrix([[1, 1,0],  [0,0,2],[0,0,1]])
print 'A = '
printmatrix(A)
print 'A^t*A = '
printmatrix(A.T*A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)

S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)
print '\n'

P = Matrix([[sqrt(2),0,0],[0,sqrt(5)/2,sqrt(5)],[0,sqrt(5)/4,sqrt(5)/2]])
print '\nP ='
printmatrix(P)

U =  Matrix([[sqrt(2)/2,sqrt(2)/2,0],[-sqrt(10)/5,sqrt(10)/5,2*sqrt(5)/5],[sqrt(10)/5,-sqrt(10)/5,sqrt(5)/5]])
print '\nU ='
printmatrix(U)

print ''

printmatrix(P*U)




print '\n\nParcial 22-06-2017. Ejercicio 5. Descomposicion en valores singulares'

A  = Matrix([[1, 1], [1,0],  [0,1]])
print 'A = '
printmatrix(A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)


S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)

print '\n\nCambio de base'
# Si T: V1 -> V2,  B1,C1 bases de V1 y B2,C2 bases de V2; entonces [T]_{C1C2} = [Id]_{C1B1} * [T]_{B1B2} * [Id]_{B2C2}
# Entonces dadas las cuatro bases y  A =  [T]_{B1B2}, calculamos    B = [T]_{C1C2}
B1 = [[4, 4, -2], [4, -2, 4], [1, -2, -2]]
C1 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
A = Matrix([[10, -2, 1], [-2, 10, 1],[-2,-2,-5]])
printmatrix(A)
B2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
C2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


U = Matrix([[0, 1, 0], [1/sqrt(2),0, 1/sqrt(2)] ,[1/sqrt(2),0, -1/sqrt(2)]])
P = Matrix([[3,0, 0], [0, 2, 0],[0,0, 1]])
printmatrix((P*U))

A =  Matrix([[2,0, 0], [0, 2, 0],[0,0, 1]])
printmatrix (U*P*U.T)


### cuentas para el final 
U = Matrix([[1/sqrt(3), -1/sqrt(3), 1/sqrt(3)], [1/sqrt(2), 1/sqrt(2),0] ,[-sqrt(6)/6, sqrt(6)/6, sqrt(6)/3]])

# print GramSchmidt([U.T*e[0], U.T*e[1], U.T*e[2]], True)

S =  Matrix([[2,0, 0], [0, 1, 0],[0,0, 0]])

V= Matrix([[0,0, 1], [-1, 0, 0],[0,1, 0]])

print ' '
printmatrix(U*S)

print ' '
printmatrix(V.T) 

print ' '
printmatrix(U*S*V.T)

A = U*S*V.T

CH = A.T*A - t*eye(3)
print ' '
printmatrix(CH)

print CH.det()

print solve(CH.det(),t)
"""

"""

A=  Matrix([[5,1, 1], [1, 5, 1],[1,1, 5]])
CH = A- t*eye(3)
print CH.det()
print teoremaVS(A)
"""

"""
# probar que T es auto adjunto
T =  Matrix([[1,1, 2], [1, 0, 1],[2, 1, 3]])
v1 = Matrix([[1], [ -1],[1]])
v2 = Matrix([[2], [ 0],[1]])
v3 = Matrix([[1], [ 1],[0]])

print T*v1
print T*v2
print T*v3


print '\n\nFinal 09-08-2017.  Descomposicion en valores singulares'

A  = Matrix([[1, 2,0],  [0,0,2],[0,0,1]])
print 'A = '
printmatrix(A)
print 'A^t*A = '
printmatrix(A.T*A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)

S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)
print '\n'


print 'Ejercicio 20.a'


print '\n\nDescomposicion polar'
A  = Matrix([[-34, 12],  [12,-41]])
print 'A = '
printmatrix(A)
print 'A^t*A = '
printmatrix(A.T*A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)

S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)
print '\n'

print 'Ejercicio 20.b'


print '\n\nDescomposicion polar'
A  = Matrix([[2, 2],  [2,-1]])
print 'A = '
printmatrix(A)
print 'A^t*A = '
printmatrix(A.T*A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)

S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)
print '\n'

print 'Ejercicio 11. matriz N'


print '\n\nDescomposicion polar'
A  = Matrix([[4, -2],  [2,-1], [0,0]])
print 'A = '
printmatrix(A)
print 'A^t*A = '
printmatrix(A.T*A)
x = teoremaVS(A)
print '\nVectores vi:',x[0],''
print 'Vectores wi:',x[1],''
print 'autovalores de A^tA:',x[2],'\n'

U = Matrix(x[1]).T
print '\nU ='
printmatrix(U)

S = Matrix(A.shape[1], A.shape[0], lambda i,j: dirac(i,j)*sqrt(x[2][i])).T
print '\nS ='
printmatrix(S)

Vt = Matrix(x[0])
print '\nV ='
printmatrix(Vt.T)

print '\nComprobando A = U*S*V^t' 
printmatrix(U*S*Vt)
print latex(U*S*Vt)

R = Rational



A= Matrix([[1,-1,2],[3,2,1],[0,1,-2]])
printmatrix(A)

print latex(A)
print det(A)
B = Matrix([[1,0,1],[R(3,4),-R(1,4),R(9,4)],[-R(3,8),R(1,8),-R(5,8)]])
printmatrix(B)
print ''
printmatrix(A.inv())

print 'aaaa'
P = Matrix([[1,1,-1],[-1,2,0],[0,1,-2]])
B = Matrix([[-1,0,0],[0,0,0],[0,0,2]])
# B = Matrix([[1,0,0],[0,2,0],[0,0,0]])
C = 5*P*B*P.inv()
print matrix(C)
print 'bbb'
print latex(C)
print det(C-t*eye(3))
print solve(det(C-t*eye(3)),t)

A = Matrix([[1,2,0],[0,0,2],[0,0,1]])
print matrix(A.T * A)
print solve(det(A.T * A-t*eye(3)),t)

A =  Matrix([[1,-3,-2,-2],[0,2,1,2],[-2,5,-1,1], [-1,3,0,1]])
A0 = A

# E(i,r,j,A):  # fila_i(A) = fila_i(A) + r*fila_j(A)
# F(r,j,A): # fila_j(A) = r*fila_j(A)
# G(i,j,A): # fila_i(A) = fila_j(A), fila_j(A) = fila_i(A)
A = E(2,2,0,eye(4))*E(3,1,0,eye(4))*A
I = E(2,2,0,eye(4))*E(3,1,0,eye(4))*eye(4)
print matrix(A) 
print matrix(I)
print '' 
A = F(-1,2,eye(4))*A
I = F(-1,2,eye(4))*I
print matrix(A) 
print matrix(I)
print '' 
A = E(0,3,2,eye(4))*E(1,-2,2,eye(4))*A
I =E(0,3,2,eye(4))*E(1,-2,2,eye(4))*I
print matrix(A) 
print matrix(I)
print ''
A = F(-R(1,2),3,eye(4))*A
I = F(-R(1,2),3,eye(4))*I
print matrix(A) 
print matrix(I)
print '' 
A = E(0,-13,3,eye(4))*E(1,9,3,eye(4))*E(2,-5,3,eye(4))*A
I = E(0,-13,3,eye(4))*E(1,9,3,eye(4))*E(2,-5,3,eye(4))*I
print matrix(A) 
print matrix(I)
print ''
A = E(0,-1,1,eye(4))*E(2,-1,1,eye(4))*E(3,-1,1,eye(4))*A
I =  E(0,-1,1,eye(4))*E(2,-1,1,eye(4))*E(3,-1,1,eye(4))*I
print matrix(A) 
print matrix(I)
print ''
A = F(2,1,eye(4))*A
I = F(2,1,eye(4))*I
print matrix(A) 
print matrix(I)
print '' 
A = G(2,3,eye(4))*G(1,2,eye(4))*A
I = G(2,3,eye(4))*G(1,2,eye(4))*I
print matrix(A) 
print matrix(I)
print '' 

print matrix(A0)
print ''
print matrix(A0*I)

print ''
print matrix(A0.inv())
print '\n'
"""

v1 = Matrix([[1,-1,2,0,1]]).T
v2 = Matrix([[2,0,1,1,0]]).T
v3 = Matrix([[0,1,0,2,0]]).T
print (3*v1-v2)
print (v1-3*v2)

print ('')
print (latex((2*v1-v3).T.tolist()))
print (latex((v1-3*v3).T.tolist()))
print (latex((v1 + v2).T.tolist()))

# E(i,r,j,A):   fila_i(A) = fila_i(A) + r*fila_j(A)
# F(r,j,A):      fila_j(A) = r*fila_j(A)
# G(i,j,A):     fila_i(A) = fila_j(A), fila_j(A) = fila_i(A)

A0 = Matrix([[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])
A = Matrix([[1,0,0,-2,-1],[0,1,0,0,-1],[0,0,1,-1,-2],[0,0,0,0,0],[0,0,0,0,0]])
A = E(0,2,1,A)
A= E(2,2,0,A)
#A= E(1,1,2,A)
#A = E(0,2,1,A)
printmatrix(A)
