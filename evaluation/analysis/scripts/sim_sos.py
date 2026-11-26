import sys, numpy as np
# usage: sim_sos.py seed out.txt   (SOS-like data: hidden LexA protein P represses 8 promoters incl. its own, G2;
# the output is the promoter activity a_i(t) (as GFP/min/OD), not P. UV cleaves P. 4 replicates, 50 points, 6 min.)
seed,out=int(sys.argv[1]),sys.argv[2]
rng=np.random.default_rng(seed)
N=8; T=50; STEP=6.0; SUB=12; dt=STEP/SUB          # Euler step 0.5 min
dP=0.02; ps=0.1                                   # protein decay (min^-1), synthesis per promoter activity
for _ in range(1000):
    tm=rng.uniform(.6,1.0,N); b=rng.uniform(.01,.04,N); K=rng.uniform(.3,1.0,N); n=rng.uniform(2,4,N)
    def act(P): return tm/(1+(P/K)**n)+b
    def sim(dose):
        P=1.0
        for _ in range(int(800/dt)): P+=dt*(ps*act(P)[1]-dP*P)           # pre-UV steady state
        X=np.zeros((N,T)); amp=0.03*np.sqrt(dose/5)*rng.uniform(.8,1.2)
        for t in range(T):
            X[:,t]=act(P)
            for s in range(SUB):
                tt=t*STEP+s*dt
                P+=dt*(ps*act(P)[1]-dP*P-amp*np.exp(-tt/60.)*P)
        return X
    runs=[sim(d) for d in (5,5,20,20)]
    allx=np.concatenate(runs,axis=1)
    if (allx.max(1)-allx.min(1)).min()>0.15: break
lo,hi=allx.min(),allx.max()
with open(out,'w') as f:
    f.write(f"repeat_number=4\ntimepoint_number={T}\n")
    for r,X in enumerate(runs):
        X=X*(1+0.03*rng.standard_normal(X.shape))
        X=np.clip(0.01+(X-lo)/(hi-lo)*0.99,0.01,1)
        for i in range(N): f.write(f"rep{r+1}\tG{i+1}\t"+"\t".join("%.6f"%v for v in X[i])+"\n")
