// Usage: node verify.js [generations]   (default 3000)
const {W,H,FIXED,initial,stepWorld,respawn,bbox,P}=require('./world.js');
const G=+(process.argv[2]||3000);
let a=initial();const hist=[a];
for(let t=1;t<=G;t++){a=stepWorld(a);respawn(a,t);hist.push(a)}
function region(g,o){const b=bbox(o[0]);let s='';for(let y=o[2]+b.y0-3;y<=o[2]+b.y1+3;y++)for(let x=o[1]+b.x0-3;x<=o[1]+b.x1+3;x++)s+=g[((y+H)%H)*W+(x+W)%W];return s}
const checks=[0.1,0.3,0.6,0.9,1].map(f=>Math.max(60,Math.floor(G*f)-5));
let fail=0;
for(const o of FIXED){const bad=checks.filter(t=>region(hist[t],o)!==region(hist[t-o[4]],o));if(bad.length)fail++;console.log(o[3].padEnd(24),bad.length?'CHANGED at '+bad:'ok')}
const band=(g,y0,y1)=>{let n=0;for(let y=y0;y<=y1;y++)for(let x=0;x<W;x++)n+=g[y*W+x];return n};
for(const[y0,y1]of[[3,40],[157,209]])console.log(`ship lanes y${y0}-${y1} pop (t / t-140):`,checks.slice(2).map(t=>band(hist[t],y0,y1)+'/'+band(hist[t-140],y0,y1)).join(' '));
// pattern self-tests
function grid(c,S,o){const g=new Uint8Array(S*S);for(const[x,y]of c)g[(y+o)*S+x+o]=1;return g}
function st(a,S){const b=new Uint8Array(S*S);for(let y=1;y<S-1;y++)for(let x=1;x<S-1;x++){let n=0;for(let j=-1;j<=1;j++)for(let i=-1;i<=1;i++)n+=a[(y+j)*S+x+i];n-=a[y*S+x];b[y*S+x]=(n===3||(n===2&&a[y*S+x]))?1:0}return b}
function period(c){const S=60,g=grid(c,S,20);let x=g;for(let p=1;p<=60;p++){x=st(x,S);if(x.every((v,i)=>v===g[i]))return p}return -1}
for(const k of['pulsar','penta','kok','queenbee','tumbler','octagon','fig8','beacon','toad','blinker'])console.log('period',k,period(P[k]));
console.log(fail?`${fail} object(s) disturbed`:'all fixed objects stable');
