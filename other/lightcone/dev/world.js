// ===== WORLD (embedded verbatim in page) =====
const W = 300, H = 210, N = W * H;
function rle(s){const out=[];let x=0,y=0,n='';for(const ch of s){if(ch>='0'&&ch<='9'){n+=ch;continue}const k=n?parseInt(n,10):1;n='';
  if(ch==='b')x+=k;else if(ch==='o'){for(let i=0;i<k;i++)out.push([x++,y])}else if(ch==='$'){y+=k;x=0}else if(ch==='!')break}return out}
function bbox(c){let a=1e9,b=1e9,d=-1e9,e=-1e9;for(const[x,y]of c){a=Math.min(a,x);b=Math.min(b,y);d=Math.max(d,x);e=Math.max(e,y)}return{x0:a,y0:b,x1:d,y1:e}}
function mirror(c,mx,my){const b=bbox(c);return c.map(([x,y])=>[mx?b.x1-x:x-b.x0,my?b.y1-y:y-b.y0])}
function transpose(c){return mirror(c.map(([x,y])=>[y,x]),false,false)}
// All orientations below were checked headlessly (period / velocity / eater absorption).
const P = {
  gun:rle('24bo$22bobo$12b2o6b2o12b2o$11bo3bo4b2o12b2o$2o8bo5bo3b2o$2o8bo3bob2o4bobo$10bo5bo7bo$11bo3bo$12b2o!'), // shoots SE
  eater:[[0,0],[1,0],[0,1],[2,1],[2,2],[2,3],[3,3]],             // absorbs the gun's gliders at gun+(90+k,76+k)
  lwss:rle('bo2bo$o4b$o3bo$4o!'),          // moves W, c/2
  mwss:rle('3bo$bo3bo$o$o4bo$5o!'),        // moves W, c/2
  hwss:rle('2b2o$o4bo$6bo$o5bo$b6o!'),     // moves E, c/2
  copperhead:rle('b2o2b2o$3b2o$3b2o$obo2bobo$o6bo2$o6bo$b2o2b2o$2b4o2$3b2o$3b2o!'), // moves N, c/10
  weekender:rle('bo12bo$bo12bo$obo10bobo$bo12bo$bo12bo$2bo3b4o3bo$6b4o$2b4o4b4o2$4bo6bo$5b2o2b2o!'), // moves N, 2c/7
  loafer:rle('b2o2bob2o$o2bo2b2o$bobo$2bo$8bo$6b3o$5bo$6bo$7b2o!'),  // moves W, c/7
  diehard:rle('6bo$2o$bo3b3o!'),           // vanishes completely after 130 gen
  pulsar:rle('2b3o3b3o2$o4bobo4bo$o4bobo4bo$o4bobo4bo$2b3o3b3o2$2b3o3b3o$o4bobo4bo$o4bobo4bo$o4bobo4bo2$2b3o3b3o!'),
  penta:rle('2bo4bo$2ob4ob2o$2bo4bo!'),
  kok:rle('2ob6o$2ob6o$2o$2o5b2o$2o5b2o$2o5b2o$7b2o$6ob2o$6ob2o!'),
  queenbee:rle('9bo$7bobo$6bobo$2o3bo2bo11b2o$2o4bobo11b2o$7bobo$9bo!'),
  tumbler:rle('bo5bo$obo3bobo$o2bobo2bo$2bo3bo$2b2ob2o!'),
  octagon:rle('3b2o$2bo2bo$bo4bo$o6bo$o6bo$bo4bo$2bo2bo$3b2o!'),
  fig8:rle('3o$3o$3o$3b3o$3b3o$3b3o!'),
  beacon:rle('2o$o$3bo$2b2o!'), toad:rle('b3o$3o!'), blinker:rle('3o!')
};
const GA={x:8,y:46}, GB={x:W-8-36,y:46}, EK=6;
// fixed objects: [cells, x, y, label, period]
const FIXED=[
  [P.gun,GA.x,GA.y,'Gosper gun',30],
  [P.eater,GA.x+90+EK,GA.y+76+EK,'eater',30],
  [mirror(P.gun,true,false),GB.x,GB.y,'Gosper gun (mirrored)',30],
  [mirror(P.eater,true,false),GB.x+35-(90+EK)-3,GB.y+76+EK,'eater',30],
  [P.queenbee,60,64,'queen bee shuttle p30',30],
  [P.pulsar,98,58,'pulsar p3',3],
  [P.tumbler,118,84,'tumbler p14',14],
  [P.fig8,118,108,'figure eight p8',8],
  [P.octagon,172,62,'octagon 2 p5',5],
  [P.penta,198,58,'pentadecathlon p15',15],
  [P.kok,176,92,"Kok's galaxy p8",8],
  [P.pulsar,20,72,'pulsar p3',3],
  [P.beacon,14,138,'beacon p2',2],
  [P.toad,34,146,'toad p2',2],
  [P.blinker,56,152,'blinker p2',2],
  [P.pulsar,268,92,'pulsar p3',3],
  [P.penta,268,128,'pentadecathlon p15',15],
  [P.kok,262,140,"Kok's galaxy p8",8]
];
// moving objects on a wrap-around (toroidal) lattice: [cells, x, y, label, vx]
const MOVERS=[
  [mirror(P.lwss,true,false),8,6,'LWSS → c/2',0.5],[mirror(P.lwss,true,false),158,6,'',0.5],
  [mirror(P.mwss,true,false),60,14,'MWSS → c/2',0.5],
  [mirror(P.hwss,true,false),W-16,23,'← HWSS c/2',-0.5],[mirror(P.hwss,true,false),130,23,'',-0.5],
  [P.lwss,W-40,33,'← LWSS c/2',-0.5],
  [mirror(transpose(P.weekender),true,false),6,160,'weekender → 2c/7',2/7],
  [transpose(P.copperhead),W-20,181,'← copperhead c/10',-0.1],
  [P.loafer,W-120,195,'← loafer c/7',-1/7]
];
// transient: diehard, optionally respawned every 300 gen
const DIEHARD=[P.diehard,12,104,'diehard (vanishes at 130)'];
function place(g,cells,ox,oy){for(const[x,y]of cells){const X=((x+ox)%W+W)%W,Y=((y+oy)%H+H)%H;g[Y*W+X]=1}}
function initial(){const g=new Uint8Array(N);for(const o of FIXED)place(g,o[0],o[1],o[2]);for(const o of MOVERS)place(g,o[0],o[1],o[2]);place(g,DIEHARD[0],DIEHARD[1],DIEHARD[2]);return g}
const XM=new Int32Array(W),XP=new Int32Array(W);for(let x=0;x<W;x++){XM[x]=(x+W-1)%W;XP[x]=(x+1)%W}
function stepWorld(a){const b=new Uint8Array(N);
  for(let y=0;y<H;y++){const up=((y+H-1)%H)*W,md=y*W,dn=((y+1)%H)*W;
    for(let x=0;x<W;x++){const l=XM[x],r=XP[x];
      const n=a[up+l]+a[up+x]+a[up+r]+a[md+l]+a[md+r]+a[dn+l]+a[dn+x]+a[dn+r];
      b[md+x]=(n===3||(n===2&&a[md+x]===1))?1:0}}
  return b}
function respawn(g,gen){if(gen%300===0)place(g,DIEHARD[0],DIEHARD[1],DIEHARD[2])}
// ===== END WORLD =====
module.exports={W,H,N,P,FIXED,MOVERS,DIEHARD,initial,stepWorld,respawn,bbox,rle,mirror,transpose};
