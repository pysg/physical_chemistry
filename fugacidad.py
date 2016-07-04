# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from numpy import linalg as LA
import pprint

#from scipy import optimize


def fugacidad(Z, P, X, R):
    '''
    Esta función que se llama fugacidad, calcula los coeficientes de fugacidad
    con una ecuación de estado para mezclas multicomponente
    T = temperatura en Kelvin
    Y = fracción molar del componente i en la fase de vapor
    X = fracción molar del componente i en la fase líquida
    '''

    T = Z[0]
    Y4 = Z[1]
    Y5 = Z[2]
    P
    Y6 = 1 - Y4 - Y5
    X4 = X[0]
    X5 = X[1]
    X6 = 1 - X4 - X5

    #print T

    #    Butano, Pentano , Hexano
    # Factor Acentrico
    w = np.array([[0.199, 0.251, 0.299]]).T
    # bar
    Pc = np.array([[38.0, 33.7, 30.1]]).T
    Pc = Pc * 100
    # K
    Tc = np.array([[425.2, 469.7, 507.5]]).T
    Tr = T / Tc
    #print "w = " , w
    #print "Pc = ", Pc
    #print "Tc = ", Tc
    #print "Tr = ", Tr
    #--------------------------------------------------------------------------
    Fw = 0.48 + (1.574 * w) - (0.176 * w ** 2)
    a = ((0.42748 * (R * Tc) ** 2) / Pc) * ((1 + Fw * (1 - (Tr ** 0.5))) ** 2)
    b = (0.08664 * R * Tc) / Pc
    #--------------------------------------------------------------------------
    #print Fw, "Parametro a:", a, b

    #Yf = np.array([Y4,Y5,Y6])
    Yf = [Y4, Y5, Y6]
    Xf = [X4, X5, X6]
    #print Yf, Xf
    #----------------------- Vapor -------------------------------------------
    amv = np.sum(Yf * a ** 0.5) ** 2
    aml = np.sum(Xf * a ** 0.5) ** 2
    #-------------------------------
    #print "amv = ", amv
    #print "aml = ", aml

    bmv = np.sum(Yf * b)
    bml = np.sum(Xf * b)
    #print "bmv = ", bmv
    #print "bml = ", bml

    Av = (amv*P)/((R*T) ** 2)
    Bv = (bmv*P)/(R*T)
    #-------------------- Liquido -------------------
    Al = (aml*P)/((R*T) ** 2)
    Bl = (bml*P)/(R*T)
    #print "Av", Av
    #print "Bv", Bv
    #print "Al", Al
    #print "Av", Bl

    Zfv = [1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]
    ZfvR = np.roots(Zfv)
    Zv = np.max(ZfvR)
    #print "Zv = ", Zv

    Zfl= [1, -1, (Al- Bl- Bl** 2), (- Al* Bl)]
    ZflR = np.roots(Zfl)
    Zl = np.min(ZflR)
    #print "Zl = ", Zl
    #------------------------------------------------------------------------------------------------------------------------------
    lnfiv = (b / bmv) * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv)  * ((b / bmv) - (2 * ((a / amv) ** 0.5))) * np.log((Zv + Bv) / Zv)
    fiv = np.exp(lnfiv)
    print "fiv = ", fiv
    lnfil = (b / bml) * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl)  * ((b / bml) - (2 * ((a / aml) ** 0.5))) * np.log((Zl + Bl) / Zl)
    fil = np.exp(lnfil)
    print "fil = ", fil
    return fil, fiv


def equilibrio(Z, P, X, R, ram):
    '''
    El metodo equilibrio calcula el equilibrio de fases liquido vapor de una mezcla multicomponente
    usando el modelo seleccionado como atributo de la clase inicial
    '''

    Fi = fugacidad(Z, P, X, R)
    #Ki = Oil / Oiv
    Ki  = Fi[0]/Fi[1]
    print "Ki = ", Ki
    Yi = np.array([Z[1], Z[2], 1 - Z[1] - Z[2]])
    Xi = np.array([X[0], X[1], 1 - X[0] - X[1]])
    print "Yi = ", Yi
    #F = Yi - Ki * Xi
    if ram == 1:
        F = Yi - np.multiply(Ki, Xi)
    elif ram == 2:
        F = Xi - np.divide(Yi, Ki)
    else:
        print "Parametro Ram no es correcto"
    print "F = ", F
    return F

def jacobianoEQ(Z, P, X, R, JEQx, ram):

    Zo = np.array([Z[0], Z[1], Z[2]])
    print "Zo", Zo
    #JEQx = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    #print "JEQx", JEQx
    dZ = 1e-3
    Zm = Zo
    for k in range(len(Zo)):
        print "k", k
        Zm[ k ] = Zo[ k ] + dZ
        Eqmas = equilibrio(Zm, P, X, R, ram)
        Zm[ k ] = Zo[ k ] - dZ
        Eqmenos = equilibrio(Zm, P, X, R, ram)
        Jc = (Eqmas - Eqmenos) / (2 * dZ)
        print "Jc[] = ", Jc[2]
        for i in range(len(Zo)):
            JEQx[i, k] = Jc[i]
    #print JEQx
    return JEQx


def newton(Z, P, X, R, ram, ep):

    Zo = Z
    print "Z = ", Z

    EqM = equilibrio(Zo, P, X, R, ram)
    print "EqM = ", EqM
    errorEq = LA.norm( EqM )

    while errorEq > ep:
        JFxo =  jacobianoEQ(Z, P, X, R, JEQx, ram)
        #print "LA.inv(JFxo) = ", LA.inv(JFxo.transpose())
        dF = LA.inv(JFxo) * 0.1 * -EqM
        dF = np.dot(LA.inv(JFxo), -EqM)
        print "dF = ", dF
        Zn = Zo + dF
        Zo = Zn
        print "Zo = ", Zo
        EqM = equilibrio(Zo, P, X, R, ram)
        errorEq = LA.norm( EqM )
    print Zo
    return EqM, JFxo, errorEq, Zn





#-------------------------------------------------
#Z = np.array([320, 0.25, 0.25]).T

Z = np.array([[320, 0.25, 0.25]]).T
#Z = Z.T

P = 400
X = np.array([[0.25, 0.25]]).T
FXo = np.array([[400, 0.25, 0.25]]).T

JEQx = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

R = 8.314
ep = 1e-4
ram = 2
#-------------------------------------------------

print "Z = ", Z
print "Z[] = ", Z[2]

Fi = fugacidad(Z, P, X, R)
print Fi

EQi = equilibrio(Z, P, X, R, ram)

print EQi

#help(fugacidad)
#help(equilibrio)

Jx = jacobianoEQ(Z, P, X, R, JEQx, ram)
#print Jx

NEqM = newton(Z, P, X, R, ram, ep)

print "Converge..."

print "\n EqM = ", NEqM[0]
print "\n JFxo = ", NEqM[1]
print "\n errorEq = ", NEqM[2]
Zn = NEqM[3]
print "\n Zn = ", Zn
print "\n T = ", Zn[0]
