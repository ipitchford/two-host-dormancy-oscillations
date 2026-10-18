from check_support import require
"""Reconstructed rigorous q-only Hopf certificate; standard library only.

Requires the exact symbolic identities certified separately by
verify_symbolic.py. No floating-point value is used in this calculation.
The exact algebraic q_H lies in the rational interval read below.
"""
import json
from fractions import Fraction as R
from pathlib import Path
from interval_arithmetic import I,CI,Dual,solve,dot,mv,DIGITS
from model import parameters,branch,jacobian,coefficients,B

LOWER='0.9769161486447415813763449883881468671926'
UPPER='0.9769161486447415813763449883881468671927'


def certificate():
    q=I.hull(LOWER,UPPER)
    a=coefficients(q);_,a1,a2,a3,a4,a5,a6=a
    T=a1*a2*a3-a3*a3+a1*a5-a1*a1*a4
    require(T.positive(), 'check failed: T.positive()')
    x=(a1*a2*a5-a3*a5-a1*a1*a6)/T
    require(x.positive(), 'check failed: x.positive()')
    omega=x.sqrt()
    quartic=[I(1),a1,a2-x,a3-a1*x,a6/x]
    _,b1,b2,b3,b4=quartic
    h1=b1;h2=b1*b2-b3;h3=b3*h2-b1*b1*b4;h4=b4*h3
    require(all(z.positive() for z in [b1,b2,b3,b4,h1,h2,h3,h4]), 'check failed: all(z.positive() for z in [b1,b2,b3,b4,h1,h2,h3,h4])')
    X=branch(q);require(all(I(z).positive() for z in X), 'check failed: all(I(z).positive() for z in X)')
    qd=Dual(q,1);pd=parameters(q=qd);Xd=branch(qd)
    Ad=jacobian(Xd,pd)
    A=[[Dual(z).value for z in row] for row in Ad]
    Ap=[[Dual(z).derivative for z in row] for row in Ad]
    iw=CI(0,omega)
    M=[[CI(A[i][j])-(iw if i==j else 0) for j in range(6)] for i in range(6)]
    pivots={}
    pivots['right_eigenvector']=[]
    u=solve([row[:5] for row in M[:5]],[-row[5] for row in M[:5]],pivots['right_eigenvector'])+[CI(1)]
    pivots['transpose_left_eigenvector']=[]
    ell=solve([[M[j][i] for j in range(5)] for i in range(5)],[-M[5][i] for i in range(5)],pivots['transpose_left_eigenvector'])+[CI(1)]
    normalisation=dot(ell,u)
    ell=[z/normalisation for z in ell]
    speed=dot(ell,mv(Ap,u)).re
    ub=[z.conj() for z in u]
    pivots['inverse_A']=[]
    h=solve(A,B(u,ub,q),pivots['inverse_A'])
    pivots['inverse_2iw_minus_A']=[]
    K=[[(2*iw if i==j else 0)-CI(A[i][j]) for j in range(6)] for i in range(6)]
    g=solve(K,B(u,u,q),pivots['inverse_2iw_minus_A'])
    normal=[-2*c+d for c,d in zip(B(u,h,q),B(ub,g,q))]
    l1=dot(ell,normal).re/(2*omega)
    require(omega.enclosed_by('0.824358242148467','0.824358242148468'), "check failed: omega.enclosed_by('0.824358242148467','0.824358242148468')")
    require(speed.enclosed_by('7.005341843234863','7.005341843234864'), "check failed: speed.enclosed_by('7.005341843234863','7.005341843234864')")
    require(l1.enclosed_by('-0.024837104528915','-0.024837104528914'), "check failed: l1.enclosed_by('-0.024837104528915','-0.024837104528914')")
    # These residual enclosures are diagnostic. Eigenvector existence follows
    # exactly from det(M)=0 and the certified invertible principal minor.
    ur=mv(M,u);lr=[dot(ell,[M[j][i] for j in range(6)]) for i in range(6)]
    require(all(z.re.contains_zero() and z.im.contains_zero() for z in ur+lr), 'check failed: all(z.re.contains_zero() and z.im.contains_zero() for z in ur+lr)')
    norm_res=dot(ell,u)-1
    require(norm_res.re.contains_zero() and norm_res.im.contains_zero(), 'check failed: norm_res.re.contains_zero() and norm_res.im.contains_zero()')
    # Independently check resident subsystem Hurwitz signs at every q in box.
    # Host 1 is independent of q. Host 2 has a=0 and b=mu3/(1-q).
    p=parameters(q=q)
    V1=R(1702,17);E1=[R(7,5),R(14,15)*V1,0,0,0,V1]
    br=R(7,5)/(1-q)
    vr=(p['l2']-1-p['C']*br)/(p['C']*br*(2*q+R(2,3)*(1-q))+1-R(1,3)*(1-q)-R(9,10)*q)
    E2=[0,0,br,2*q*vr*br,R(2,3)*(1-q)*vr*br,vr]
    def char(mat):
        n=len(mat);W=[[I(int(i==j)) for j in range(n)] for i in range(n)];out=[I(1)]
        for k in range(1,n+1):
            AW=[[sum((mat[i][h]*W[h][j] for h in range(n)),I(0)) for j in range(n)] for i in range(n)]
            c=-sum((AW[i][i] for i in range(n)),I(0))/k;out.append(c)
            W=[[AW[i][j]+(c if i==j else 0) for j in range(n)] for i in range(n)]
        return out
    resident=[]
    for E,ids in [(E1,[0,1,5]),(E2,[2,3,4,5])]:
        J=jacobian(E,p);r=char([[I(J[i][j]) for j in ids] for i in ids])
        ds=[r[1],r[1]*r[2]-r[3]]
        if len(r)==5:ds.extend([r[3]*ds[1]-r[1]**2*r[4],r[4]*(r[3]*ds[1]-r[1]**2*r[4])])
        else:ds.append(r[3]*ds[1])
        require(all(z.positive() for z in r+ds), 'check failed: all(z.positive() for z in r+ds)')
        resident.append({'coefficients':[z.to_json() for z in r],'Hurwitz_determinants':[z.to_json() for z in ds]})
    # Effective invader active growth after eliminating infected/dormant states.
    H1=sum(E1[:5]);H2=sum(E2[:5])
    invasion2=p['l2']-1-p['C']*H1-V1+R(1,3)*(1-q)*V1+R(9,10)*q*V1
    invasion1=p['l1']-1-p['C']*H2-R(2,3)*vr
    require(invasion1.positive() and invasion2.positive(), 'check failed: invasion1.positive() and invasion2.positive()')
    return {'status':'PASS','provenance':'fresh reconstruction, 17 October 2026; no lost-file byte-identity claim',
      'decimal_grid_digits':DIGITS,'q_interval':[LOWER,UPPER],
      'omega':omega.to_json(),'total_real_crossing_speed':speed.to_json(),'first_Lyapunov_coefficient':l1.to_json(),
      'x_omega_squared':x.to_json(),'D3':T.to_json(),
      'stable_quartic_coefficients':[z.to_json() for z in quartic],
      'stable_quartic_Hurwitz_determinants':[z.to_json() for z in [h1,h2,h3,h4]],
      'equilibrium':[I(z).to_json() for z in X],
      'right_eigenvector':[z.to_json() for z in u],
      'transpose_left_eigenvector':[z.to_json() for z in ell],
      'raw_left_right_pairing':normalisation.to_json(),
      'certified_elimination_pivots':pivots,
      'resident_stability':resident,
      'resident2_viral_density':vr.to_json(),
      'effective_invasion_growth':[invasion1.to_json(),invasion2.to_json()],
      'normalisation':'u_6=1; ell^T u=1; A u=i omega u; ell^T A=i omega ell^T',
      'scope':'local supercritical Hopf with stable positive cycles for q>q_H sufficiently close; no certified radius'}


if __name__=='__main__':
    result=certificate()
    dest=Path(__file__).with_name('hopf_interval_certificate.json')
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','decimal_grid_digits','omega','total_real_crossing_speed','first_Lyapunov_coefficient']},indent=2))
