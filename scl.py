#rref or ref
n2 = int(input("Enter number of equations:"))
n1 = int(input("Enter number of unknowns:"))
augmatrix = []
for i in range(n2):
    print("Equation "+str(i+1)+": ")
    a = []
    for j in range(n1):
        n = int(input("Enter coefficient of x"+str(j+1)+": "))
        a.append(n)
    b = int(input("Enter RHS of the eqn: "))
    a.append(b)
    augmatrix.append(a)
print("The augmented matrix is:")
for i in range(n2):
    print(augmatrix[i])

pos = -1;
flag = 0;
zero = [0]*(n1+1);
for i in range(n2):
    if augmatrix[i] == zero:
        for k in range(i+1,n2):
            if augmatrix[k] != zero:
                print("Not a REF")
                flag = 1
    else:
        for index in range(0,n1+1):
            if augmatrix[i][index] != 0:
                if index > pos:
                    pos = index                    
                    break;
                else:
                    print("Neither REF nor RREF");
                    flag = 1
                    break;

                    
        for k in range(i+1,n2):
            if augmatrix[k][pos] != 0:
                print("Neither REF nor RREF");
                flag = 1;
                break;
    
    if flag == 1:
        break;
        
        
    else:
        if augmatrix[i][pos] != 1:
            print(pos)
            print('The matrix is in REF but not in RREF')
            break;
        
        for k in range(0,n2):
            if k!=i:
                print(pos)
                if augmatrix[k][pos] !=0:
                    flag = 1;
                    print("The matrix is in REF but not in RREF");
                    break;
         
    if flag == 1:
        break;




#Gauss jacobi

from pprint import pprint
from numpy import array, zeros, diag, diagflat, dot

def jacobi(A, b, N=25, x=None):
    if x is None:
        x = zeros(len(A[0]))
    D = diag(A)
    R = A - diagflat(D)
    for i in range(N):
        x = (b - dot(R, x)) / D
    return x


n = int(input("Enter the number of equations: "))
m = int(input("Enter the number of unknowns: "))
a = []
c=[]
for i in range(n):
    b = []
    for j in range(m+1):
        if j == m:
            num = int(input("Enter constant term: "))
            c.append(num)
        else:
            num = int(input("Enter matrix entry ["+str(i)+"]["+str(j)+"]: "))
            b.append(num)
    a.append(b)
# print(a)
# print(c)
print("Input Matrix :")
for i in range(len(a)):
    print(a[i])
    
A = array(a)
b = array(c)
guess = array([1.0]*n)
sol = jacobi(A, b, N=25, x=guess)
print("A: ")
pprint(A)
print("B: ")
pprint(b)
print("x = ")
pprint(sol)


#Gauss jordan elimination
n = int(input("Enter the number of rows of matrix: "))
m = int(input("Enter the number of cols of matrix: "))

a = []
final = []
for i in range(n):
    b = []
    for j in range(m):
        num = int(input("Enter matrix entry ["+str(i)+"]["+str(j)+"]: "))
        b.append(num)
    a.append(b);
    
    
print("Input Matrix :")
for i in range(len(a)):
    print(a[i])

flag = 1;
iternum = 0
zero = [0]*m
b=[]
while flag:
    c = 0
    b = []
    for i in a:
        if i != zero:
            b.append(i)
        else:
            c +=1
    j=0
    while j<c:
        b.append(zero)
        j+=1
    
    pos = m+1;
    value = -1000
    swap = 0;
    for i in range(n-c):
        for j in range(m):
            if b[i][j] != 0:
                if j<=pos and value<b[i][j]:
                    value = b[i][j]
                    pos = j
                    swap = i
                    
            break;
    if pos == m+1:
        pos = 0
    
    c_ = b[swap]
    b[swap] = b[0]
    b[0] = c_
    
    #print(b[0])
    
    scalar = b[0][pos]
    
    if scalar!=0:
        for i in range(m):
            b[0][i] /=scalar
    
    print("After Iteration " + str(iternum) + " : ")
    for i in range(len(b)):
        print(b[i])
    
    for i in range(1,n):
        scalar = b[i][pos]
        if scalar != 0:
            for j in range(m):
                b[i][j] -= scalar*b[0][j]
     
    #print("After making all the values in the column of leading 1 zero : ",b)
    
    if len(b) == 1:
        flag = 0
        d = [0]*iternum + b[0]
        final.append(d)
        
    else:
        d = [0]*iternum + b[0]
        final.append(d)
        a = []
        for i in range(1,n):
            a_ = []
            for j in range(1,m):
                a_.append(b[i][j])
            a.append(a_)
        
        #print("Matrix for next iteration : ",a)
        n = n-1
        m = m-1
        iternum +=1 
        
        
        
print("Final matrix after performing Gauss Jordan Elimination : ")
for i in range(len(final)):
    print(final[i])


#Gauss seidel
import re 

def isDiagonalDominant(a):
    for i in range(len(a)):
        sum = 0
        for j in range(len(a[0])):
            if i!= j:
                sum += a[i][j]
        if a[i][i] >= sum:
            continue
        else:
            return False
    return True

n = int(input("Enter the number of equations: "))
eqns = []
for i in range(n):
    s = input("Enter equation " + str(i+1) + ":")
    eqns.append(s)
a = []
for x in eqns:
    x = re.findall("(-?[\d]+)", x)
    x = [int(i) for i in x]
    a.append(x)
    
print("Matrix form: ")
for i in range(len(a)):
    print(a[i])


if(isDiagonalDominant(a)):
    i = 0
    initail  = [0]*n
    values = [0]*n
    maxiter = 25
    while i < maxiter:
        for j in range(len(a)):
            k =0
            values[j] = a[j][-1]
            while k < (len(a[j])-1):
                if k!=j:
                    values[j] -= values[k]*a[j][k]
                k+=1      
            values[j] /= a[j][j];
        print(values)
        i+=1;
        

#diagonalize
import sympy as sp
import numpy as np
x=sp.var('x')

def eigen_vector(mat,l):
    length=mat.shape[0]
    I=sp.eye(length)
    Z=I.multiply(l)
    H=mat-Z
    J=H.nullspace()
    return J
    
def eigen_value(mat):
    length=mat.shape[0]
    I=sp.eye(length)
    Z=I.multiply(x)
    H=mat-Z
    Det=H.det()
    print("Diagonlisation of {} is \n".format(mat.tolist()))
    L=sp.solve(Det,x)
    l=[]
    for i in L:
        k=eigen_vector(mat,i)

        for j in range(len(k)):
            l.append(k[j].tolist()) 
    M=np.empty((mat.shape[0],len(l)))
    for i in range(len(l)):
        for j in range(mat.shape[0]):
            M[j][i]=l[i][j][0]
    
    Mi=np.linalg.inv(M)
    print(M)
    K=Mi.dot(mat)
    K=K.dot(M)
    print(K)
       
A = sp.Matrix([[1,0],[6,-1]])
eigenvalue(A)
p,d = A.diagonalize()
print(p)


#eigen valiues
from sympy import var,Matrix,shape,eye,solve,multiply,det

x = var('x')


A = Matrix([[4,0,1],[-1,2,0],[2,0,1]])
l = shape(A)
I = eye(l[0])
H = A-(I*x)
det = det()
eigenvalues = solve(det,x)
print("Eigen values of A are ",eigenvalues)
for e in eigenvalues:
    H = A - (I.multiply(e))
    J=H.nullspace()
    for i in range(len(J)):
        print("Eigen vector for eigen value ",e," is: ",J[i].tolist())


#Inverse power method

import numpy as np


def inversePower(a,x):
    for i in range(10):
        x = np.dot(a, x)
        x_n = x / x.max()
    print('Inverse Power method:')
    print('Eigenvector of matrix')
    for i in range(len(a)):
        print(a[i])
    print('is ',x_n)
    find_ray(a,x_n)
    
def find_ray(a,x_n):
    h=np.dot(a,x_n)
    m=np.dot(h,x_n)
    
    x_2=np.dot(x_n,x_n)
    dom_eival=m/x_2
    print('Reyleigh quotient method:')
    print('The smallest eigen value of matrix {} is {:.2f}'.format(a.tolist(),dom_eival))
    

a=np.array([[2,-1,0],[-1,2,-1],[0,-1,2]])  
x=np.array([1,1,1]) 
inversePower(a,x)


#power method
import numpy as np


def power(a,x):
    for i in range(10):
        x = np.dot(a, x)
        x_n = x / x.min()
    print('Power method:')
    print('Eigentvector of matrix')
    for i in range(len(a)):
        print(a[i])
    print("is ", x_n)
    find_ray(a,x_n)
    
def find_ray(a,x_n):
    h=np.dot(a,x_n)
    m=np.dot(h,x_n)
    
    x_2=np.dot(x_n,x_n)
    dom_eival=m/x_2
    print('Reyleigh quotient method:')
    print('The dominant eigen value for matrix {} is {:.2f}'.format(a.tolist(),dom_eival))    


ab = np.array([[2, -12], 
              [1, -5]])   
x=np.array([1,1]) 
power(ab,x)


#Bisection
import math

def f(x):
    return 3*x + math.cos(x) - x
  
def bisection(a,b):
 
    if (f(a) * f(b) >= 0):
        print("You have not assumed right a and b\n")
        return
  
    c = a
    while ((b-a) >= 0.01):
        c = (a+b)/2
        if (f(c) == 0.0):
            break
        if (f(c)*f(a) < 0):
            b = c
        else:
            a = c
             
    print("The value of root is : ","%.4f"%c)
     
# Driver code
# Initial values assumed
a =-1
b = 0
bisection(a, b)

#Fixed point
from sympy import symbols,sympify,diff

x = symbols('x')

def func(equation,k):
    return equation.subs(x,k).evalf()

def diff_eval(equation,k):
    return func(diff(equation,x),k)

def fixed_point(equation,a,b,f):
    if func(equation,a)*func(equation,b)>0:
        print("Give a correct value for a and b")
    else:
        difa = abs(diff_eval(f,a));
        difb = abs(diff_eval(f,b));

        print(difa,"\n");
        print(difb,"\n");

        if max(difa,difb)>=1:
            print("Convergence criteria not satisfied");
        else:
            while abs(func(f,a)-a) > 0.0001:
                a = func(f,a)
            print("The root is :", "%.4f"%a)

def main():
    equation = input("Enter the equation")
    equation = sympify(equation)

    print(equation)

    a = int(input("Enter A : "))
    b = int(input("Enter B : "))

    f = sympify(input("Enter an iteration function"));

    fixed_point(equation,a,b,f);

main()


#newton raphson

def f(x):
    return x * x * x - x * x + 2

def df(x):
    return 3 * x * x - 2 * x

def newtonRaphson( x ):
    h = f(x) / df(x)
    while abs(h) >= 0.0001:
        h = f(x)/df(x)
        # x(i+1) = x(i) - f(x) / f'(x)
        x = x - h
    print("The value of the root is : ", "%.4f"% x)
 
x0 = 0.1 # Initial value
newtonRaphson(x0)

#regula falsi
from sympy import symbols, sympify
x = symbols('x')

def func(val,k): 
    return k.subs(x,val).evalf()

def regulaFalsi(k,a,b):
 
    if (func(a,k) * func(b,k) >= 0):
        print("Invalid values of a and b")
        return
  
    c = a 
    for i in range(1000):
        c = (a*func(b,k)-b*func(a,k))/(func(b,k)-func(a,k))
        if (func(c,k) == 0.0):
            break
        if (func(c,k)*func(a,k) < 0):
            b = c
        else:
            a = c
             
    print("The value of root is : ","%.4f"%c)
     
k = input("Enter the equation : ")
k = sympify(k)
print(k)

a = int(input("Enter the value of a : "))
b = int(input("Enter the value of b : "))
regulaFalsi(k,a,b)
