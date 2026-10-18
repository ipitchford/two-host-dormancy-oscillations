// Run with Node.js. t = 2*C - 1 is nonnegative on the parameter ray.
const certificate = (() => {
const add=(a,b)=>Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||0n)+(b[i]||0n));
const mul=(a,b)=>{const c=Array(a.length+b.length-1).fill(0n);a.forEach((x,i)=>b.forEach((y,j)=>c[i+j]+=x*y));return c;};
const neg=a=>a.map(x=>-x);
const det=A=>{
  let out=[0n];
  const rec=(row,used,product,cols)=>{
    if(row===A.length){
      let inv=0;
      for(let i=0;i<cols.length;i++)for(let j=i+1;j<cols.length;j++)inv+=cols[i]>cols[j];
      out=add(out,inv%2?neg(product):product);return;
    }
    for(let col=0;col<A.length;col++)if(!used.has(col)&&A[row][col].some(x=>x!==0n)){
      used.add(col);rec(row+1,used,mul(product,A[row][col]),[...cols,col]);used.delete(col);
    }
  };
  rec(0,new Set(),[1n],[]);
  while(out.length>1&&out.at(-1)===0n)out.pop();
  return out;
};
// Coefficients of 36000*p_((t+1)/2)(s), descending in s.
// Each array lists polynomial coefficients in t in increasing order.
const a=[[36000n],[824760n,378000n],[6517148n,4951380n],
[19960218n,18467484n],[20883247n,19832547n],
[13223868n,11876130n],[12122649n,12122649n]];
const result=[];
for(let k=1;k<=6;k++){
  const D=det(Array.from({length:k},(_,i)=>Array.from({length:k},(_,j)=>a[2*j-i+1]||[0n])));
  const desired=k<=4?1n:-1n;
  if(!D.every(x=>desired*x>0n))throw Error("Hurwitz sign certificate failed at k="+k);
  result.push({k,degree:D.length-1,sign:desired===1n?"+":"-",coefficients:D.map(String)});
}
return {status:"PASS",hurwitz:result};
})();
console.log(JSON.stringify(certificate,null,2));
