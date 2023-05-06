from math import factorial, sqrt, isclose

def wigner3j(j1, j2, j3, m1, m2, m3):
    '''Returns the square of the specified Wigner 3-j symbol'''

    CG = 0
    for k in range(max(0, j3 - m3 + m1 - j1), min(j2 - m2, j3 - m3) + 1): # range gotten by examining inequalities obtained from the six factorial factors in summand
        CG += ((-1)**(j1-m1-j3+m3+k)*(factorial(j1+m1+j3-m3-k)*factorial(j2+m2+k)) / 
               (factorial(k)*factorial(j3-m3-k)*factorial(j1-m1-j3+m3+k)*factorial(j2-m2-k)))

    CG **= 2
    
    CG *= ((factorial(j1+j2-j3)*(2*j3+1)*factorial(j3+m3)*factorial(j3-m3)*factorial(j1-m1)*factorial(j2-m2)) / 
           (factorial(j1-j2+j3)*factorial(j2-j1+j3)*factorial(j1+j2+j3+1)*factorial(j1+m1)*factorial(j2+m2)))
    
    return 1/(2*j3+1)*CG
