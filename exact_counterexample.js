const certificate = (() => {
// Exact rational certificate for the dormancy counterexample.
// Run with Node.js. No packages, numerical eigenvalue solver or file I/O are used.
const gcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b];}return a;};
const R=(n,d=1n)=>{n=BigInt(n);d=BigInt(d);if(d===0n)throw Error("division by zero");if(d<0n){n=-n;d=-d;}const g=gcd(n,d);return [n/g,d/g];};
const Z=R(0), O=R(1);
const add=(a,b)=>R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const neg=a=>[-a[0],a[1]], sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>R(a[0]*b[0],a[1]*b[1]), div=(a,b)=>R(a[0]*b[1],a[1]*b[0]);
const sum=a=>a.reduce(add,Z), eq=(a,b)=>a[0]===b[0]&&a[1]===b[1];
const show=a=>a[1]===1n?String(a[0]):a[0]+"/"+a[1];
const sign=a=>a[0]>0n?"+":a[0]<0n?"-":"0";
const assert=(p,m)=>{if(!p)throw Error(m);};
const mm=(A,B)=>A.map(row=>B[0].map((_,j)=>sum(row.map((x,k)=>mul(x,B[k][j])))));
const characteristic=A=>{
  const n=A.length; let B=A.map((r,i)=>r.map((_,j)=>i===j?O:Z));const out=[O];
  for(let k=1;k<=n;k++){const AB=mm(A,B),c=neg(div(sum(AB.map((r,i)=>r[i])),R(k)));out.push(c);B=AB.map((r,i)=>r.map((x,j)=>i===j?add(x,c):x));}return out;
};
const determinant=A=>{
  const B=A.map(r=>r.slice()); let d=O;
  for(let k=0;k<B.length;k++){
    let p=k;while(p<B.length&&B[p][k][0]===0n)p++;
    if(p===B.length)return Z;
    if(p!==k){[B[p],B[k]]=[B[k],B[p]];d=neg(d);}
    const pivot=B[k][k];d=mul(d,pivot);
    for(let i=k+1;i<B.length;i++){const t=div(B[i][k],pivot);for(let j=k+1;j<B.length;j++)B[i][j]=sub(B[i][j],mul(t,B[k][j]));B[i][k]=Z;}
  }return d;
};
const hurwitz=a=>Array.from({length:a.length-1},(_,k)=>determinant(Array.from({length:k+1},(_,i)=>Array.from({length:k+1},(_,j)=>a[2*j-i+1]||Z))));
const routh=a=>{const n=a.length-1,m=Math.ceil((n+1)/2);let prev=Array.from({length:m},(_,j)=>a[2*j]||Z),cur=Array.from({length:m},(_,j)=>a[2*j+1]||Z),first=[prev[0],cur[0]];
for(let k=2;k<=n;k++){assert(cur[0][0]!==0n,"zero Routh pivot");const next=Array.from({length:m},(_,j)=>div(sub(mul(cur[0],prev[j+1]||Z),mul(prev[0],cur[j+1]||Z)),cur[0]));first.push(next[0]);prev=cur;cur=next;}return first;};
const parameters=C=>({C,mu:R(1),D:R(1),v:R(1),r:R(1,2),q:R(49,50),k:R(1,20),sigma:R(9,20),m:R(3),mu3:R(7,5),l1:div(add(R(13),mul(R(665),C)),R(3)),l2:add(R(467,300),div(mul(R(665),C),R(3)))});
const field=(X,p)=>{
const [a,i,b,d,j,V]=X,H=sum(X.slice(0,5)),T=mul(p.C,H),DV=mul(p.D,V),u=sub(O,p.q),rv=add(p.r,p.v),ks=add(mul(p.k,p.mu),p.sigma);
return [add(mul(a,sub(sub(sub(p.l1,p.mu),T),DV)),mul(p.r,i)),sub(mul(DV,a),mul(rv,i)),
sum([mul(b,sub(sub(sub(p.l2,p.mu),T),DV)),mul(p.r,j),mul(p.sigma,d)]),
sub(mul(mul(p.q,DV),b),mul(ks,d)),sub(mul(mul(u,DV),b),mul(rv,j)),
sub(sub(mul(mul(p.m,p.v),add(i,j)),mul(DV,add(a,mul(u,b)))),mul(p.mu3,V))];};
const jacobian=(X,p)=>{
const [a,i,b,d,j,V]=X,H=sum(X.slice(0,5)),T=mul(p.C,H),DV=mul(p.D,V),Ca=mul(p.C,a),Cb=mul(p.C,b),u=sub(O,p.q),rv=add(p.r,p.v),ks=add(mul(p.k,p.mu),p.sigma);
return [
[sub(sub(sub(sub(p.l1,p.mu),T),DV),Ca),sub(p.r,Ca),neg(Ca),neg(Ca),neg(Ca),neg(mul(p.D,a))],
[DV,neg(rv),Z,Z,Z,mul(p.D,a)],
[neg(Cb),neg(Cb),sub(sub(sub(sub(p.l2,p.mu),T),DV),Cb),sub(p.sigma,Cb),sub(p.r,Cb),neg(mul(p.D,b))],
[Z,Z,mul(p.q,DV),neg(ks),Z,mul(mul(p.q,p.D),b)],
[Z,Z,mul(u,DV),Z,neg(rv),mul(mul(u,p.D),b)],
[neg(DV),mul(p.m,p.v),neg(mul(u,DV)),Z,mul(p.m,p.v),neg(add(mul(p.D,add(a,mul(u,b))),p.mu3))]
];};
const X=[R(1),R(10,3),R(20),R(196),R(4,3),R(5)];
const expected=C=>[[18000,0],[223380,378000],[782884,4951380],[746367,18467484],[525350,19832547],[673869,11876130],[0,12122649]].map(([a,b])=>div(add(R(a),mul(R(b),C)),R(18000)));
const interpolation=[];
for(let j=1;j<=7;j++){const C=R(j,2),p=parameters(C);assert(field(X,p).every(x=>x[0]===0n),"nonzero equilibrium residual");const a=characteristic(jacobian(X,p));assert(a.every((x,i)=>eq(x,expected(C)[i])),"characteristic mismatch");interpolation.push(show(C));}
// Each characteristic coefficient has degree at most six in C,
// since every Jacobian entry is affine in C. Seven points prove equality.
const p=parameters(R(1,2)),V1=R(1702,17),V2=R(38195,34589);
const E1=[R(7,5),mul(R(14,15),V1),Z,Z,Z,V1];
const E2=[Z,Z,R(70),mul(R(686,5),V2),mul(R(14,15),V2),V2];
const choose=(A,ix)=>ix.map(i=>ix.map(j=>A[i][j]));
const boundary=[];
for(const [name,E,ix] of [["host1",E1,[0,1,5]],["host2",E2,[2,3,4,5]]]){
  assert(field(E,p).every(x=>x[0]===0n),name+" residual");
  const ar=characteristic(choose(jacobian(E,p),ix)),ds=hurwitz(ar);
  assert(ds.every(x=>x[0]>0n),name+" resident stability");
  const af=characteristic(jacobian(E,p)),first=routh(af);
  const signs=first.map(sign),changes=signs.slice(1).reduce((n,s,i)=>n+(s!==signs[i]),0);
  assert(!signs.includes("0")&&changes===1,name+" invasion index");
  boundary.push({name,residentCoefficients:ar.map(show),residentHurwitzDeterminants:ds.map(show),fullRouthSigns:signs});
}
const a=characteristic(jacobian(X,p)),hf=hurwitz(a),rs=routh(a).map(sign);
assert(hf.map(sign).join("")==="++++--","coexistence determinant signs");
assert(rs.join("")==="+++++-+","coexistence Routh signs");
assert(eq(sum(jacobian(X,p).map((r,i)=>r[i])),R(-2291,100)),"trace check");
return {status:"PASS",exactInterpolationPoints:interpolation,coexistenceJacobian:jacobian(X,p).map(r=>r.map(show)),coexistencePolynomialTimes36000:a.map(x=>show(mul(x,R(36000)))),coexistenceHurwitzSigns:hf.map(sign),coexistenceRouthSigns:rs,boundary,trace:"-2291/100",H:"665/3"};
})();
console.log(JSON.stringify(certificate,null,2));
