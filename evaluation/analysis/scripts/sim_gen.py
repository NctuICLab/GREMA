import sys, numpy as np
# usage: gen.py seed reps noise out.txt   (IRMA topology, HFODE, Euler delta_t=1 like GREMA)
seed,reps,noise,out=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),sys.argv[4]
T=16; N=5
edges={(1,2):1,(2,4):1,(2,3):1,(3,5):1,(3,1):1,(3,4):1,(4,2):-1,(5,1):-1}
rng=np.random.default_rng(seed)
while True:
    tm=rng.uniform(.4,1,N); b=rng.uniform(0,.05,N); deg=rng.uniform(.2,.5,N)
    k=rng.uniform(.2,.6,(N,N)); n=rng.uniform(2,4,(N,N))
    def sim(x0):
        X=np.zeros((N,T)); X[:,0]=x0
        for t in range(T-1):
            for i in range(N):
                p=1.0
                for j in range(N):
                    if (j+1,i+1) in edges:
                        e=n[i,j]*edges[(j+1,i+1)]
                        p*=1/(1+(k[i,j]/X[j,t])**e)
                dx=tm[i]*p+b[i]-deg[i]*X[i,t]
                X[i,t+1]=min(max(X[i,t]+dx,0.01),1.5)
        return X
    runs=[sim(rng.uniform(.05,.9,N)) for _ in range(reps)]
    allx=np.concatenate(runs,axis=1)
    if (allx.max(1)-allx.min(1)).min()>0.25: break   # every gene has visible dynamics
lo,hi=allx.min(),allx.max()
with open(out,'w') as f:
    f.write(f"repeat_number={reps}\ntimepoint_number={T}\n")
    for r,X in enumerate(runs):
        X=X*(1+noise*rng.standard_normal(X.shape))
        X=0.01+(X-lo)/(hi-lo)*0.99; X=np.clip(X,0.01,1)
        for i in range(N): f.write(f"rep{r+1}\tG{i+1}\t"+"\t".join("%.6f"%v for v in X[i])+"\n")
