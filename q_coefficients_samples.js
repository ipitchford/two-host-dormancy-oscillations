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

const polynomial=(a,q)=>a.reduceRight((v,c)=>add(mul(v,q),R(c)),Z);
const power=(a,n)=>Array.from({length:n}).reduce(v=>mul(v,a),O);
const coefficients=q=>{
 const L=add(mul(R(5),q),R(49)),q2=power(q,2),L2=power(L,2);
 return [O,
 div(polynomial([48020,140287,184865],q),mul(R(300),mul(q,L))),
 div(polynomial([1176490,7087752,32467841,47850335],q),mul(R(9000),mul(q2,L))),
 div(polynomial([-39883011,1273247899,31184629,1429336425],q),mul(R(90000),mul(q2,L))),
 div(polynomial([269383385929,-3765832488152,29089900503650,-22001310797200,607919428125],q),mul(R(2700000),mul(q2,L2))),
 div(mul(R(-49),polynomial([-1292021318,21642947739,-25203074295,-19568133750,14410621250],q)),mul(R(450000),mul(q2,L2))),
 div(mul(R(833),mul(polynomial([-833,17020],q),polynomial([833,-16824,17125],q))),mul(R(600000),mul(q2,L)))
 ];
};
const checks=[];
for(let k=0;k<=12;k++){
 const q=R(970+k,1000),p=parameters(R(1,2));p.q=q;
 const V=div(R(49),mul(R(10),q));
 const b=div(sub(mul(R(17020),q),R(833)),mul(R(15),mul(q,add(mul(R(5),q),R(49)))));
 const a=sub(R(7,5),mul(sub(O,q),b)),i=mul(R(2,3),mul(V,a)),d=mul(R(2),mul(q,mul(V,b))),j=mul(R(2,3),mul(sub(O,q),mul(V,b)));
 const X=[a,i,b,d,j,V],actual=characteristic(jacobian(X,p)),expected=coefficients(q);
 assert(field(X,p).every(z=>z[0]===0n),"q-branch residual");
 assert(actual.every((z,i)=>eq(z,expected[i])),"q coefficient mismatch");
 checks.push(show(q));
}
return {status:"PASS",q_values:checks,scope:"Exact branch residuals and all six characteristic coefficients at thirteen rational q values."};
})();
console.log(JSON.stringify(certificate,null,2));
