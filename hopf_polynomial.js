const certificate = (() => {
const trim=a=>{a=a.slice();while(a.length>1&&a.at(-1)===0n)a.pop();return a;};
const add=(a,b)=>trim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||0n)+(b[i]||0n)));
const mul=(a,b)=>{const c=Array(a.length+b.length-1).fill(0n);a.forEach((x,i)=>b.forEach((y,j)=>c[i+j]+=x*y));return trim(c);};
const scale=(a,k)=>a.map(x=>x*k),neg=a=>scale(a,-1n);
const gcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b];}return a;};
const det=A=>{let out=[0n];const rec=(row,used,product,cols)=>{if(row===A.length){let inv=0;for(let i=0;i<cols.length;i++)for(let j=i+1;j<cols.length;j++)inv+=cols[i]>cols[j];out=add(out,inv%2?neg(product):product);return;}for(let col=0;col<A.length;col++)if(!used.has(col)&&A[row][col].some(x=>x!==0n)){used.add(col);rec(row+1,used,mul(product,A[row][col]),[...cols,col]);used.delete(col);}};rec(0,new Set(),[1n],[]);return trim(out);};
const ev=(a,n,d)=>{let out=a.at(-1),dp=d;for(let k=a.length-2;k>=0;k--){out=out*n+a[k]*dp;dp*=d;}return out;};
const divideLinear=a=>{const r=a.slice(),b=Array(a.length-1).fill(0n);for(let i=a.length-1;i>=1;i--){if(r[i]%5n!==0n)throw Error("nonintegral quotient");b[i-1]=r[i]/5n;r[i]-=5n*b[i-1];r[i-1]-=49n*b[i-1];}if(r[0]!==0n)throw Error("nonzero remainder");return trim(b);};
const q=[0n,1n],L=[49n,5n],W=scale(mul(mul(q,q),mul(L,L)),5400000n);
const a=[W,
scale(mul(mul(q,L),[48020n,140287n,184865n]),18000n),
scale(mul(L,[1176490n,7087752n,32467841n,47850335n]),600n),
scale(mul(L,[-39883011n,1273247899n,31184629n,1429336425n]),60n),
scale([269383385929n,-3765832488152n,29089900503650n,-22001310797200n,607919428125n],2n),
scale([-1292021318n,21642947739n,-25203074295n,-19568133750n,14410621250n],-588n),
scale(mul(L,mul([-833n,17020n],[833n,-16824n,17125n])),7497n)];
let N=det(Array.from({length:5},(_,i)=>Array.from({length:5},(_,j)=>a[2*j-i+1]||[0n])));
const rawDegree=N.length-1;let qFactors=0,LFactors=0;
while(N[0]===0n){N=N.slice(1);qFactors++;}
while(ev(N,-49n,5n)===0n){N=divideLinear(N);LFactors++;}
const content=N.reduce((g,x)=>gcd(g,x),0n),sgn=N.at(-1)>0n?1n:-1n;
const P=N.map(x=>x/content*sgn);
if(P.length!==18)throw Error("Unexpected primitive numerator degree");
const lower=9769161486447415813763449883881468671926n,upper=lower+1n,den=10n**40n;
const lo=ev(P,lower,den),hi=ev(P,upper,den);
if(lo*hi>=0n)throw Error("No root sign change");
const D=P.slice(1).map((x,i)=>BigInt(i+1)*x),midN=lower+upper,midD=2n*den;
const derivativeNumerator=ev(D,midN,midD);
const M2=P.reduce((s,x,k)=>s+BigInt(k*(k-1))*(x<0n?-x:x),0n);
const derivativeMargin=(derivativeNumerator<0n?-derivativeNumerator:derivativeNumerator)*2n*den-M2*(midD**16n);
if(derivativeMargin<=0n)throw Error("Derivative sign not isolated");
return {status:"PASS",rawNumeratorDegree:rawDegree,cancelled_q_factors:qFactors,cancelled_5q_plus_49_factors:LFactors,primitiveDegree:17,coefficients_in_ascending_order:P.map(String),rootLower:"0.9769161486447415813763449883881468671926",rootUpper:"0.9769161486447415813763449883881468671927",P_at_lower_sign:lo>0n?"+":"-",P_at_upper_sign:hi>0n?"+":"-",derivative_sign_on_interval:derivativeNumerator>0n?"+":"-",uniqueness_verified:true};
})();
console.log(JSON.stringify(certificate,null,2));
