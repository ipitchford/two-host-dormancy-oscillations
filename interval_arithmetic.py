from check_support import require
"""Transparent outward integer intervals on a fixed decimal grid.

Only Python's unbounded integers, Fraction, and math.isqrt affect bounds.
Every arithmetic operation encloses the exact set operation. Complex
division excludes zero using the enclosure of squared modulus. Midpoints
are used only to order candidate pivots; all accepted pivots are certified.
"""
from fractions import Fraction
from math import isqrt

DIGITS=60
SCALE=10**DIGITS


def ceildiv(a,b):
    return -((-a)//b)


class I:
    __slots__=('lo','hi')
    def __init__(self, value=0):
        if isinstance(value,I): self.lo,self.hi=value.lo,value.hi; return
        f=Fraction(value)
        self.lo=(f.numerator*SCALE)//f.denominator
        self.hi=ceildiv(f.numerator*SCALE,f.denominator)
    @classmethod
    def raw(cls,lo,hi):
        require(lo<=hi, 'check failed: lo<=hi')
        z=cls.__new__(cls);z.lo,z.hi=lo,hi;return z
    @classmethod
    def hull(cls,lo,hi):
        return cls.raw(I(lo).lo,I(hi).hi)
    def __add__(self,other):
        if isinstance(other,(CI,Dual)):return NotImplemented
        b=I(other);return I.raw(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I.raw(-self.hi,-self.lo)
    def __sub__(self,other):return self+-I(other) if not isinstance(other,(CI,Dual)) else NotImplemented
    def __rsub__(self,other):return I(other)+-self
    def __mul__(self,other):
        if isinstance(other,(CI,Dual)):return NotImplemented
        b=I(other);v=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return I.raw(min(v)//SCALE,ceildiv(max(v),SCALE))
    __rmul__=__mul__
    def __truediv__(self,other):
        if isinstance(other,(CI,Dual)):return NotImplemented
        b=I(other)
        if b.lo<=0<=b.hi:raise ArithmeticError('real interval denominator contains zero')
        pairs=[(a*SCALE,c) for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return I.raw(min(a//c for a,c in pairs),max(ceildiv(a,c) for a,c in pairs))
    def __rtruediv__(self,other):return I(other)/self
    def square(self):
        high=max(self.lo*self.lo,self.hi*self.hi)
        low=0 if self.lo<=0<=self.hi else min(self.lo*self.lo,self.hi*self.hi)
        return I.raw(low//SCALE,ceildiv(high,SCALE))
    def __pow__(self,n):
        require(isinstance(n,int) and n>=0, 'check failed: isinstance(n,int) and n>=0')
        if n==0:return I(1)
        if n==1:return self
        if n%2:return self*self**(n-1)
        return (self**(n//2)).square()
    def sqrt(self):
        if self.lo<0:raise ArithmeticError('negative square-root interval')
        lo=isqrt(self.lo*SCALE);hi=isqrt(self.hi*SCALE)
        if hi*hi<self.hi*SCALE:hi+=1
        return I.raw(lo,hi)
    def positive(self):return self.lo>0
    def negative(self):return self.hi<0
    def contains_zero(self):return self.lo<=0<=self.hi
    def enclosed_by(self,lower,upper):
        return Fraction(self.lo,SCALE)>Fraction(lower) and Fraction(self.hi,SCALE)<Fraction(upper)
    def to_json(self):
        def fmt(n):
            s='-' if n<0 else '';n=abs(n)
            return s+str(n//SCALE)+'.'+str(n%SCALE).zfill(DIGITS)
        return [fmt(self.lo),fmt(self.hi)]
    def __repr__(self):return '['+', '.join(self.to_json())+']'


class CI:
    __slots__=('re','im')
    def __init__(self,re=0,im=0):
        if isinstance(re,CI):self.re,self.im=re.re,re.im;return
        self.re,self.im=I(re),I(im)
    def __add__(self,b):b=CI(b);return CI(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return CI(-self.re,-self.im)
    def __sub__(self,b):return self+-CI(b)
    def __rsub__(self,b):return CI(b)+-self
    def __mul__(self,b):
        b=CI(b);return CI(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=CI(b);den=b.re.square()+b.im.square()
        if not den.positive():raise ArithmeticError('complex interval denominator not separated from zero')
        return CI((self.re*b.re+self.im*b.im)/den,(self.im*b.re-self.re*b.im)/den)
    def __rtruediv__(self,b):return CI(b)/self
    def conj(self):return CI(self.re,-self.im)
    def midpoint_score(self):return (self.re.lo+self.re.hi)**2+(self.im.lo+self.im.hi)**2
    def to_json(self):return {'real':self.re.to_json(),'imaginary':self.im.to_json()}
    def __repr__(self):return str(self.to_json())


class Dual:
    """Forward automatic differentiation with outward interval components."""
    __slots__=('value','derivative')
    def __init__(self,value,derivative=0):
        if isinstance(value,Dual):self.value,self.derivative=value.value,value.derivative;return
        self.value,self.derivative=I(value),I(derivative)
    def __add__(self,b):b=Dual(b);return Dual(self.value+b.value,self.derivative+b.derivative)
    __radd__=__add__
    def __neg__(self):return Dual(-self.value,-self.derivative)
    def __sub__(self,b):return self+-Dual(b)
    def __rsub__(self,b):return Dual(b)+-self
    def __mul__(self,b):
        b=Dual(b);return Dual(self.value*b.value,self.derivative*b.value+self.value*b.derivative)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=Dual(b);return Dual(self.value/b.value,(self.derivative*b.value-self.value*b.derivative)/(b.value*b.value))
    def __rtruediv__(self,b):return Dual(b)/self
    def __pow__(self,n):
        require(isinstance(n,int) and n>=0, 'check failed: isinstance(n,int) and n>=0')
        out=Dual(1)
        for _ in range(n):out=out*self
        return out


def solve(A,b,pivot_log=None):
    n=len(A);M=[[CI(x) for x in row]+[CI(v)] for row,v in zip(A,b)]
    require(all(len(row)==n+1 for row in M), 'check failed: all(len(row)==n+1 for row in M)')
    for k in range(n):
        candidates=sorted(range(k,n),key=lambda i:M[i][k].midpoint_score(),reverse=True)
        p=next((i for i in candidates if (M[i][k].re.square()+M[i][k].im.square()).positive()),None)
        if p is None:raise ArithmeticError('no certified nonzero pivot')
        M[k],M[p]=M[p],M[k]
        pivot=M[k][k]
        if pivot_log is not None:pivot_log.append(pivot.to_json())
        for i in range(k+1,n):
            f=M[i][k]/pivot
            for j in range(k+1,n+1):M[i][j]=M[i][j]-f*M[k][j]
            M[i][k]=CI(0)
    x=[CI(0) for _ in range(n)]
    for i in range(n-1,-1,-1):
        x[i]=(M[i][n]-sum((M[i][j]*x[j] for j in range(i+1,n)),CI(0)))/M[i][i]
    return x


def dot(x,y):return sum((a*b for a,b in zip(x,y)),CI(0))
def mv(A,x):return [dot(row,x) for row in A]


if __name__=='__main__':
    import json
    # Cross-check directed arithmetic against exact rational corner operations.
    vals=[Fraction(-7,3),Fraction(-1,7),Fraction(0),Fraction(2,11),Fraction(19,5)]
    cases=0
    for a in vals:
        for b in vals:
            for op,exact in [(lambda x,y:x+y,a+b),(lambda x,y:x-y,a-b),(lambda x,y:x*y,a*b)]:
                v=op(I(a),I(b));require(Fraction(v.lo,SCALE)<=exact<=Fraction(v.hi,SCALE), 'check failed: Fraction(v.lo,SCALE)<=exact<=Fraction(v.hi,SCALE)');cases+=1
            if b:
                v=I(a)/I(b);exact=a/b
                require(Fraction(v.lo,SCALE)<=exact<=Fraction(v.hi,SCALE), 'check failed: Fraction(v.lo,SCALE)<=exact<=Fraction(v.hi,SCALE)');cases+=1
    for a in [Fraction(0),Fraction(1,7),Fraction(2),Fraction(1001,3)]:
        v=I(a).sqrt();require(Fraction(v.lo,SCALE)**2<=a<=Fraction(v.hi,SCALE)**2, 'check failed: Fraction(v.lo,SCALE)**2<=a<=Fraction(v.hi,SCALE)**2');cases+=1
    print(json.dumps({'status':'PASS','exact_rational_corner_tests':cases,'decimal_grid_digits':DIGITS}))
