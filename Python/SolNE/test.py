from solne import sne_ud_1, sne_ud_2, sne_fd_1

f = ''
x0 = 0
tol = 0.0001

#------------------------------------------------------------------------------
# Test 1: sne_ud_1 (Euler's Method)
#------------------------------------------------------------------------------
f = '(x**2) - exp(x) - 3*x + 2'; x0 = 0.7
sne_ud_1(f, x0, tol)

#------------------------------------------------------------------------------
# Test 2:
#------------------------------------------------------------------------------
