# -----------------------------------------------------------------------------------------------------
# General complexities: C_alg = 2^(k*w+d)
# -----------------------------------------------------------------------------------------------------    
def gb_compl_kd(dreg, nv, ne):
    # Assumption: nv=ne
    # Returns: k,d such that C_GB = 2^(k*w+d)
    
    if dreg == None:
        return None, None
    
    # Formula used for C_GB: binomial(nv + dreg, nv)**w
    k = float(log(binomial(ne + dreg, ne),2))
    d = 0
    
    k = round(ceil(k * 10) / 10, 1)
    d = round(ceil(d * 10) / 10, 1)

    return k,d

def fglm_compl_kd(dI, nv, ne):
    # Assumption: nv=ne
    # Returns: k,d such that C_FGLM = 2^(k*w+d)
    
    if dI == None:
        return None, None
    
    # Formula used for C_FGLM: nv * dI**w.
    k = float(log(dI,2))
    d = float(log(nv,2))
    
    k = round(ceil(k * 10) / 10, 1)
    d = round(ceil(d * 10) / 10, 1)
    
    return k,d

def print_stats(eqs,w_gb=2.37,w_fglm=3):
    print('='*80)
    eqs = Sequence(eqs)
    neqs = len(eqs)
    nvars = Sequence(eqs).nvariables()
    print(f"{neqs} equations in {nvars} variables", flush=True)
    
    print('-'*80)
    mac = macaulay_bound(eqs)
    k,d = gb_compl_kd(mac, nvars, neqs)
    print(f"Macaulay bound:  {mac}", flush=True)
    print(f"GB complexity:   {k:.1f}w + {d:.1f} ({round(k*w_gb+d,1)})", flush=True)
    
    print('-'*80)
    bez = bezout_bound(eqs)
    k,d = fglm_compl_kd(bez, nvars, neqs)
    print(f"Bézout bound:    {factor(bez)}", flush=True)
    print(f"FGLM complexity: {k:.1f}w + {d:.1f} ({round(k*w_fglm+d,1)})", flush=True)
    print('='*80)