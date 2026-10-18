"""One explicit transcription of the six-dimensional dormancy vector field.

Reconstructed from the durable research Page on 17 October 2026.
Arithmetic is generic: Fraction, SymPy expressions, or outward intervals.
State order is (a, i, b, d, j, V); q is the dormancy probability.
"""
from fractions import Fraction as R


def parameters(C=R(1, 2), q=R(49, 50)):
    return dict(C=C, q=q, mu=1, D=1, v=1, r=R(1, 2),
                k=R(1, 20), sigma=R(9, 20), m=3, mu3=R(7, 5),
                l1=(13+665*C)/3, l2=R(467, 300)+665*C/3)


def field(X, p):
    a, i, b, d, j, V = X
    H = sum(X[:5]); DV = p['D']*V; C = p['C']; q = p['q']
    r = p['r']; v = p['v']; mu = p['mu']; sig = p['sigma']
    return [a*(p['l1']-mu-C*H-DV)+r*i,
            DV*a-(r+v)*i,
            b*(p['l2']-mu-C*H-DV)+r*j+sig*d,
            q*DV*b-(p['k']*mu+sig)*d,
            (1-q)*DV*b-(r+v)*j,
            -DV*(a+(1-q)*b)+p['m']*v*(i+j)-p['mu3']*V]


def jacobian(X, p):
    a, i, b, d, j, V = X
    H = sum(X[:5]); C = p['C']; D = p['D']; q = p['q']
    r = p['r']; v = p['v']; mu = p['mu']; sig = p['sigma']
    Ca=C*a; Cb=C*b; DV=D*V
    return [
      [p['l1']-mu-C*H-DV-Ca, r-Ca, -Ca, -Ca, -Ca, -D*a],
      [DV, -(r+v), 0, 0, 0, D*a],
      [-Cb, -Cb, p['l2']-mu-C*H-DV-Cb, sig-Cb, r-Cb, -D*b],
      [0, 0, q*DV, -(p['k']*mu+sig), 0, q*D*b],
      [0, 0, (1-q)*DV, 0, -(r+v), (1-q)*D*b],
      [-DV, p['m']*v, -(1-q)*DV, 0, p['m']*v,
       -D*(a+(1-q)*b)-p['mu3']]]


def branch(q):
    V=49/(10*q)
    b=(17020*q-833)/(15*q*(5*q+49))
    a=R(7,5)-(1-q)*b
    return [a, R(2,3)*V*a, b, 2*q*V*b, R(2,3)*(1-q)*V*b, V]


def coefficients(q):
    """Recovered monic characteristic coefficients, descending in z."""
    L=5*q+49
    return [1,
      (184865*q**2+140287*q+48020)/(300*q*L),
      (47850335*q**3+32467841*q**2+7087752*q+1176490)/(9000*q**2*L),
      (1429336425*q**3+31184629*q**2+1273247899*q-39883011)/(90000*q**2*L),
      (607919428125*q**4-22001310797200*q**3+29089900503650*q**2
       -3765832488152*q+269383385929)/(2700000*q**2*L**2),
      -49*(14410621250*q**4-19568133750*q**3-25203074295*q**2
           +21642947739*q-1292021318)/(450000*q**2*L**2),
      833*(17020*q-833)*(17125*q**2-16824*q+833)/(600000*q**2*L)]


def B(x,y,q,C=R(1,2)):
    """The symmetric bilinear Hessian D^2 F, with D=1."""
    hxy=x[0]*sum(y[:5])+y[0]*sum(x[:5])
    kxy=x[2]*sum(y[:5])+y[2]*sum(x[:5])
    av=x[0]*y[5]+y[0]*x[5]
    bv=x[2]*y[5]+y[2]*x[5]
    return [-C*hxy-av, av, -C*kxy-bv, q*bv, (1-q)*bv,
            -av-(1-q)*bv]
