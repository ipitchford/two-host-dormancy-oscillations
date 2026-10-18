from check_support import require
"""Fresh exact symbolic verification of the recovered dormancy results.

SymPy 1.14.0 is used for rational identities only. The interval normal-form
calculation is separate and uses Python's standard library exclusively.
"""
import json,time
from pathlib import Path
import sympy as s
from model import parameters,field,jacobian,branch,coefficients,B

HERE=Path(__file__).resolve().parent


def cleaned_char(A,z):
    # Newton/Faddeev recurrence with rational cancellation at each step.
    A=s.Matrix(A);n=A.rows;W=s.eye(n);out=[s.Integer(1)]
    for k in range(1,n+1):
        AW=(A*W).applyfunc(s.cancel)
        a=s.cancel(-s.trace(AW)/k);out.append(a)
        W=AW+a*s.eye(n)
    return out


def hurwitz(a,k):
    def entry(i,j):
        h=2*j-i+1
        return a[h] if 0<=h<len(a) else 0
    return s.cancel(s.Matrix(k,k,entry).det(method='domain-ge'))


def positive_ray(expr,C,t,desired=1):
    num,den=s.fraction(s.cancel(expr.subs(C,(t+1)/2)))
    pn,pd=s.Poly(num,t),s.Poly(den,t)
    if pd.LC()<0:pn,pd=-pn,-pd
    require(all(c>0 for c in pd.all_coeffs()), 'denominator coefficient sign')
    require(all(desired*c>0 for c in pn.all_coeffs()), 'numerator coefficient sign')
    return {'numerator_coefficients_descending':[str(c) for c in pn.all_coeffs()],
            'denominator_coefficients_descending':[str(c) for c in pd.all_coeffs()],
            'strict_sign':'positive' if desired>0 else 'negative'}


def main():
    start=time.monotonic();C,t,q,z=s.symbols('C t q z')
    p=parameters(C=C);X=s.Matrix([1,s.Rational(10,3),20,196,s.Rational(4,3),5])
    symbols=s.symbols('a i b d j V');F=s.Matrix(field(symbols,p))
    symbolic_J=F.jacobian(symbols)
    require(symbolic_J==s.Matrix(jacobian(symbols,p)), 'check failed: symbolic_J==s.Matrix(jacobian(symbols,p))')
    require(all(s.cancel(v)==0 for v in field(X,p)), 'check failed: all(s.cancel(v)==0 for v in field(X,p))')
    J=s.Matrix(jacobian(X,p));a=cleaned_char(J,z)
    expected=[1,(223380+378000*C)/s.Integer(18000),(782884+4951380*C)/s.Integer(18000),
              (746367+18467484*C)/s.Integer(18000),(525350+19832547*C)/s.Integer(18000),
              (673869+11876130*C)/s.Integer(18000),12122649*C/s.Integer(18000)]
    require(all(s.cancel(v-w)==0 for v,w in zip(a,expected)), 'check failed: all(s.cancel(v-w)==0 for v,w in zip(a,expected))')
    family={'positive_equilibrium_residuals_zero':True,'jacobian_is_exact_derivative':True,
            'characteristic_coefficients':[str(v) for v in a],
            'Hurwitz_sign_certificates':[],'resident_certificates':[]}
    for k in range(1,7):
        d=hurwitz(a,k);entry=positive_ray(d,C,t,1 if k<=4 else -1);entry['k']=k
        family['Hurwitz_sign_certificates'].append(entry)
    V1=(25+1652*C)/(5+7*C);V2=5*(167+45500*C)/(167+207200*C)
    E1=[s.Rational(7,5),s.Rational(14,15)*V1,0,0,0,V1]
    E2=[0,0,70,s.Rational(686,5)*V2,s.Rational(14,15)*V2,V2]
    for name,E,ids,other in [('host1',E1,[0,1,5],[2,3,4]),('host2',E2,[2,3,4,5],[0,1])]:
        require(all(s.cancel(v)==0 for v in field(E,p)), 'check failed: all(s.cancel(v)==0 for v in field(E,p))')
        EJ=s.Matrix(jacobian(E,p));r=cleaned_char(EJ.extract(ids,ids),z)
        result={'name':name,'resident_state':[str(v) for v in E],
                'resident_characteristic_coefficients':[str(v) for v in r],
                'resident_Hurwitz_sign_certificates':[],
                'invader_matrix':[[str(EJ[i,j]) for j in other] for i in other]}
        for k in range(1,len(r)):
            entry=positive_ray(hurwitz(r,k),C,t);entry['k']=k
            result['resident_Hurwitz_sign_certificates'].append(entry)
        family['resident_certificates'].append(result)
    family['invasion_sign_certificates']=[positive_ray(V1-5,C,t),positive_ray(5-V2,C,t)]
    # Virus-free host thresholds determine the remaining unstable dimensions.
    K1=(p['l1']-1)/C;K2=(p['l2']-1)/C
    family['virus_free_infection_thresholds']=[positive_ray(K1-s.Rational(7,5),C,t),positive_ray(K2-70,C,t)]
    (HERE/'family_symbolic_certificate.json').write_text(json.dumps(family,indent=2)+'\n')
    print('all-C family and resident certificates: PASS',flush=True)

    pq=parameters(q=q);Xq=branch(q)
    require(all(s.cancel(v)==0 for v in field(Xq,pq)), 'check failed: all(s.cancel(v)==0 for v in field(Xq,pq))')
    Aq=s.Matrix(jacobian(Xq,pq)).applyfunc(s.cancel)
    aq=cleaned_char(Aq,z);target=coefficients(q)
    require(all(s.cancel(v-w)==0 for v,w in zip(aq,target)), 'check failed: all(s.cancel(v-w)==0 for v,w in zip(aq,target))')
    a0,a1,a2,a3,a4,a5,a6=aq
    T=s.cancel(a1*a2*a3-a3*a3+a1*a5-a1*a1*a4)
    omega2=s.cancel((a1*a2*a5-a3*a5-a1*a1*a6)/T)
    D5=hurwitz(aq,5)
    numerator,denominator=s.fraction(D5)
    _,P=s.Poly(numerator,q).primitive()
    if P.LC()<0:P=-P
    recovered=[int(v) for v in json.loads((HERE/'hopf_polynomial_coefficients.json').read_text())]
    require(list(reversed([int(v) for v in P.all_coeffs()]))==recovered, 'check failed: list(reversed([int(v) for v in P.all_coeffs()]))==recovered')
    require(P.degree()==17, 'check failed: P.degree()==17')
    odd=s.cancel(a1*omega2**2-a3*omega2+a5)
    even=s.cancel(-omega2**3+a2*omega2**2-a4*omega2+a6)
    quotients={}
    for name,expr in [('odd',odd),('even',even)]:
        num,den=s.fraction(expr);poly=s.Poly(num,q)
        quotient,remainder=poly.div(P)
        require(remainder.is_zero, 'check failed: remainder.is_zero')
        quotients[name]={'numerator_quotient_by_P':str(quotient.as_expr()),'denominator':str(den)}
    # Every possible denominator factor is q, 5q+49, or numerator of T.
    # All are excluded from zero by the outward certificate at the isolated root.
    beta=[s.Integer(1),a1,s.cancel(a2-omega2),s.cancel(a3-a1*omega2),s.cancel(a6/omega2)]
    product=s.Poly((z*z+omega2)*sum(beta[k]*z**(4-k) for k in range(5)),z)
    factor_identity=[]
    for k,(lhs,rhs) in enumerate(zip(product.all_coeffs(),aq)):
        num,den=s.fraction(s.cancel(lhs-rhs))
        if num==0:factor_identity.append({'coefficient':k,'identically_zero':True})
        else:
            quotient,remainder=s.Poly(num,q).div(P);require(remainder.is_zero, 'check failed: remainder.is_zero')
            factor_identity.append({'coefficient':k,'zero_modulo_P':True,'denominator':str(den)})
    # Symbolic Hessian comparison ensures the normal form's B convention.
    xs=s.symbols('x0:6');ys=s.symbols('y0:6');pqF=s.Matrix(field(symbols,pq))
    computedB=[]
    for fk in pqF:
        Hess=s.hessian(fk,symbols)
        computedB.append((s.Matrix(xs).T*Hess*s.Matrix(ys))[0])
    require(all(s.expand(v-w)==0 for v,w in zip(computedB,B(xs,ys,q))), 'check failed: all(s.expand(v-w)==0 for v,w in zip(computedB,B(xs,ys,q)))')
    require(all(s.diff(f,symbols[i],symbols[j],symbols[k])==0 for f in pqF for i in range(6) for j in range(6) for k in range(6)), 'check failed: all(s.diff(f,symbols[i],symbols[j],symbols[k])==0 for f in pqF for i in range(6) for j in range(6) for k in range(6))')
    hopf={'status':'PASS','branch_residuals_identically_zero':True,
          'characteristic_identity_for_all_q':True,'characteristic_coefficients':[str(v) for v in aq],
          'primitive_degree':17,'primitive_P_coefficients_ascending':[str(v) for v in recovered],
          'D5_denominator_factored':str(s.factor(denominator)),
          'D5_numerator_multiplier_of_P':str(s.cancel(numerator/P.as_expr())),
          'omega_squared':str(s.factor(omega2)),'D3':str(s.factor(T)),
          'even_odd_remainders':quotients,'quadratic_quartic_factorisation':factor_identity,
          'quartic_coefficients':[str(s.factor(v)) for v in beta],
          'B_is_exact_Hessian':True,'third_derivative_identically_zero':True,
          'jacobian_branch':[[str(v) for v in Aq.row(i)] for i in range(6)],
          'jacobian_total_q_derivative':[[str(s.factor(s.diff(v,q))) for v in Aq.row(i)] for i in range(6)]}
    (HERE/'hopf_symbolic_certificate.json').write_text(json.dumps(hopf,indent=2)+'\n')
    print(json.dumps({'status':'PASS','family':'all C >= 1/2','q_branch':'rational-function identity',
                      'Hopf_P_degree':17,'quadratic_factor':'exact modulo P','runtime_seconds':time.monotonic()-start},indent=2),flush=True)


if __name__=='__main__':main()
