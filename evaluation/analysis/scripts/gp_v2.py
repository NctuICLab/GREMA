import sys, numpy as np
# GP augmentation v2. Hyperparameters are marginalised over a grid (weights = marginal likelihood).
# usage: gpv2.py input M outdir seed noisy(0/1). No gold standard is used.
inp,M,outdir,seed,noisy=sys.argv[1],int(sys.argv[2]),sys.argv[3],int(sys.argv[4]),int(sys.argv[5])
L=open(inp).read().splitlines(); T=int(L[1].split('=')[1]); real=L[2:]
t=np.arange(T,dtype=float); D=np.abs(t[:,None]-t[None,:])
def kern(name,l):
    r=D/l
    if name=="rbf": return np.exp(-0.5*r**2)
    if name=="m52": return (1+np.sqrt(5)*r+5*r**2/3)*np.exp(-np.sqrt(5)*r)
    return (1+np.sqrt(3)*r)*np.exp(-np.sqrt(3)*r)
names=("rbf","m52","m32"); ls=np.geomspace(1,T,14); nzs=np.geomspace(1e-3,1,12); s2s=np.geomspace(0.3,3,5)
rng=np.random.default_rng(seed)
samples=[]; info=[]
for ln in real:
    p=ln.split('\t'); y0=np.array(p[2:],float); mu0,sd0=y0.mean(),y0.std(); y=(y0-mu0)/sd0
    grid=[];lml=[]
    for nm in names:
        for l in ls:
            Kf=kern(nm,l)
            for s2 in s2s:
                for nz in nzs:
                    K=s2*Kf+nz*np.eye(T)
                    sign,ld=np.linalg.slogdet(K)
                    lml.append(-0.5*y@np.linalg.solve(K,y)-0.5*ld-0.5*T*np.log(2*np.pi)); grid.append((nm,l,s2,nz))
    lml=np.array(lml); w=np.exp(lml-lml.max()); w/=w.sum()
    ess=1/np.sum(w**2); info.append((p[1],ess))
    out=[]
    for m in range(M):
        nm,l,s2,nz=grid[rng.choice(len(grid),p=w)]
        Kf=s2*kern(nm,l); K=Kf+nz*np.eye(T); Ki=np.linalg.inv(K)
        mu=Kf@Ki@y; cov=Kf-Kf@Ki@Kf; cov=(cov+cov.T)/2+1e-9*np.eye(T)
        f=rng.multivariate_normal(mu,cov)
        if noisy: f=f+np.sqrt(nz)*rng.standard_normal(T)
        out.append(np.clip(f*sd0+mu0,0.01,1.0))
    samples.append(np.array(out))
for m in range(M):
    with open(f"{outdir}/in_{m}.txt","w") as f:
        f.write(f"repeat_number=1\ntimepoint_number={T}\n")
        for i,ln in enumerate(real): f.write(f"rep1\t{info[i][0]}\t"+"\t".join("%.6f"%v for v in samples[i][m])+"\n")
print(" ".join(f"{g}:ESS={e:.0f}" for g,e in info))
