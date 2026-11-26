import sys, numpy as np
# Fully data-driven GP: kernel type, length scale and noise are chosen by closed-form leave-one-out
# log predictive density. Search ranges come from the data (length scale: 1 step .. T steps).
# No gold standard is used anywhere in this file.
inp,M,outdir=sys.argv[1],int(sys.argv[2]),sys.argv[3]
L=open(inp).read().splitlines(); T=int(L[1].split('=')[1]); real=L[2:]
t=np.arange(T,dtype=float); D=np.abs(t[:,None]-t[None,:])
def kern(name,l):
    r=D/l
    if name=="rbf": return np.exp(-0.5*r**2)
    if name=="m52": return (1+np.sqrt(5)*r+5*r**2/3)*np.exp(-np.sqrt(5)*r)
    if name=="m32": return (1+np.sqrt(3)*r)*np.exp(-np.sqrt(3)*r)
def loo(K,y):
    Ki=np.linalg.inv(K); a=Ki@y; d=np.diag(Ki)
    mu=y-a/d; var=1/d
    lpd=-0.5*np.log(2*np.pi*var)-0.5*(y-mu)**2/var
    cover=np.mean(np.abs(y-mu)<=1.96*np.sqrt(var))
    return lpd.sum(),cover
rng=np.random.default_rng(int(sys.argv[4]))
samples=[];info=[]
for ln in real:
    p=ln.split('\t'); y0=np.array(p[2:],float); mu0,sd0=y0.mean(),y0.std(); y=(y0-mu0)/sd0
    best=None
    for name in ("rbf","m52","m32"):
        for l in np.geomspace(1,T,14):
            for nz in np.geomspace(1e-3,1,12):
                K=kern(name,l)*(1-0)+nz*np.eye(T)   # unit signal variance (y is standardised)
                s,c=loo(K,y)
                if best is None or s>best[0]: best=(s,name,l,nz,c)
    s,name,l,nz,c=best; info.append((p[1],name,l,nz,c))
    Kf=kern(name,l); K=Kf+nz*np.eye(T); Ki=np.linalg.inv(K)
    post_mu=Kf@Ki@y; post_cov=Kf-Kf@Ki@Kf
    post_cov=(post_cov+post_cov.T)/2+1e-9*np.eye(T)
    f=rng.multivariate_normal(post_mu,post_cov,size=M)+np.sqrt(nz)*rng.standard_normal((M,T))
    samples.append(np.clip(f*sd0+mu0,0.01,1.0))
for m in range(M):
    with open(f"{outdir}/in_{m}.txt","w") as f:
        f.write(f"repeat_number=1\ntimepoint_number={T}\n")
        for i,ln in enumerate(real): f.write(f"rep1\t{info[i][0]}\t"+"\t".join("%.6f"%v for v in samples[i][m])+"\n")
for g,name,l,nz,c in info: print(f"{g} kernel={name} length={l:.2f} noise={nz:.3f} LOO95%coverage={c:.2f}")
