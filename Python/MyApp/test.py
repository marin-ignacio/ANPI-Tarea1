from solne import sne_ud_1, sne_ud_2,sne_ud_3,sne_ud_4,sne_ud_5,sne_ud_6, sne_fd_1, sne_fd_2, sne_fd_3,sne_fd_4,sne_fd_5,sne_fd_6

f = ''
x0 = 0
tol = 0.0001
ans = []

#------------------------------------------------------------------------------
# TEST OF DERIVATIVE METHODS
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Test 1: sne_ud_1 (Euler's Method)
#------------------------------------------------------------------------------
f = '(x**2) - exp(x) - 3*x + 2'; x0 = 0.7

#Expected Answer: x=0.25752586
ans = sne_ud_1(f, x0, tol)

print('root: x=0.25752586')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 2: sne_ud_2 (Dong's method)
#------------------------------------------------------------------------------
f = '(exp((x**2) + (7*x) - 30) - 1)**6'; x0 = 3.1; m = 6

#Expected Answer: x=3.00722843
ans = sne_ud_2(f, m, x0, tol)

print('root: x=3.00722843')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 3: sne_ud_3 (Homeier's method)
#------------------------------------------------------------------------------
f = 'x**2'; x0 = 1; m = 5

#Expected Answer: x = 0.038709845246447605
ans = sne_ud_3(f, m, x0, tol,graf)

print('root: x = 0.025780961748129938')
print('x_aprox:', ans[0])
print('iter:', ans[9])
print()

#------------------------------------------------------------------------------
# Test 4: sne_ud_4 (Ostrowski)
#------------------------------------------------------------------------------
f = 'x**2-13*x-12'; x0 = 2; tol = 0.01

#Expected Answer: x = -0.8656255778919385
ans = sne_ud_4(f, x0, tol)

print('root: x=-0.8656255778919385')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 5: sne_ud_5 (Halley)
#------------------------------------------------------------------------------
f = 'e**x-(1/x)'; x0 = 1; tol = 0.001

#Expected Answer: x = 0.5673043631887621
ans = sne_ud_4(f, x0, tol)

print('root: x=0.5673043631887621')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 6: sne_ud_6 (Chebyshev)
#------------------------------------------------------------------------------
f = 'e**x-(1/x)'; x0 = 1; tol = 0.001

#Expected Answer: x = 0.5671397801525238
ans = sne_ud_6(f, x0, tol)

print('root: x=0.5671397801525238')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# TEST OF FREE DERIVATIVE METHODS
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Test 1: sne_fd_1 (Improved Otrowski's Method Free Derivative)
#------------------------------------------------------------------------------
f = '(x**3) + 4*(x**2) - 10'; x0 = 1.5

#Expected Answer: x=1.36523001
ans = sne_fd_1(f, x0, tol)

print('root: x=1.36523001')
print('x_aprox:', ans[0])
print('iter:', ans[1])

#------------------------------------------------------------------------------
# Test 2: sne_fd_2 (Jain method based on Steffensen secant method)
#------------------------------------------------------------------------------
f = '(x**2+5*x'; x0 = 1;tol = 0.001

#Expected Answer: x=0.08235294117647074
ans = sne_fd_2(f, x0, tol,graf)

print('root: x=3.00722843')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 3: sne_fd_3 (Dehghan and Hajarian method)
#------------------------------------------------------------------------------
f = 'cos(x)'; x0 = 1;tol = 0.001

#Expected Answer: x = 0.5403023058681398
ans = sne_fd_3(f, x0, tol,graf)

print('root: x=3.00722843')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 4: sne_fd_4 (Ostrowski method)
#------------------------------------------------------------------------------
f = 'x**2-13*x-12'; x0 = 2;tol = 0.01

#Expected Answer: x = -0.86567552952761451
ans = sne_fd_4(f, x0, tol)

print('root: x=-0.86567552952761451')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()

#------------------------------------------------------------------------------
# Test 5: sne_fd_5 (Muller Bisection method)
#------------------------------------------------------------------------------
f = 'x**2-13*x-12'; x0 = 3; x1 = 1; x3 = 2

#Expected Answer: x = 2.618033988749895
ans = sne_fd_5(f, x0, x1, x2)

print('root: x=2.618033988749895')
print('x_aprox:', ans[0])
print()

#------------------------------------------------------------------------------
# Test 6: sne_fd_6 (Steffensen)
#------------------------------------------------------------------------------
f = 'e**x-(1/x)'; x0 = 1; tol = 0.001

#Expected Answer: x = 0.5671399085157954
ans = sne_fd_6(f, x0, tol)

print('root: x=0.5671399085157954')
print('x_aprox:', ans[0])
print('iter:', ans[1])
print()
