#__        __    _ _         __  __    _    ____ __  __    _       __ _ _
#\ \      / / __(_) |_ ___  |  \/  |  / \  / ___|  \/  |  / \     / _(_) | ___
# \ \ /\ / / '__| | __/ _ \ | |\/| | / _ \| |  _| |\/| | / _ \   | |_| | |/ _ \
#  \ V  V /| |  | | ||  __/ | |  | |/ ___ \ |_| | |  | |/ ___ \  |  _| | |  __/
#   \_/\_/ |_|  |_|\__\___| |_|  |_/_/   \_\____|_|  |_/_/   \_\ |_| |_|_|\___|


# ==============================================================================
# General part
# ==============================================================================
def gen_magma_file(filename, eqs, fixed_elems, mdir):
    """
    Generates a Magma file for Gröbner basis computation on a given system of equations.

    Parameters:
    - filename (str): The name of the Magma file to be created.
    - eqs (list): A list of polynomial equations to be used in the computation.
    - fixed_elems (list): Variables that have previously been fixed; these are added as comment in the file.
    - mdir (str): The directory path where the generated file will be saved.
    
    Returns:
    None
      The function writes the equations and configurations needed for Magma's Gröbner basis and variety computations to the specified file.
    """
    F = eqs[0].base_ring()
    p = F.characteristic()
    P = eqs[0].parent()
    
    with open(f"{mdir}/{filename}", 'w') as magma_file:
        # Write fixed elements as comment
        magma_file.write("// {};\n".format(fixed_elems))
        magma_file.write("// {};\n".format([f.degree() for f in eqs]))
        
        print("Writing eqns to " + magma_file.name)
        magma_file.write("p := {};\n".format(p))
        magma_file.write("R<{}> := PolynomialRing(GF(p),{}, \"grevlex\");\n".format(",".join(str(i) for i in P.gens()), P.ngens()))
        magma_file.write("I := ideal<R |\n")
        for i, f in enumerate(eqs):
            if i < len(eqs) - 1:
                magma_file.write("{},\n".format(f))
            else:
                magma_file.write("{}\n".format(f))
        magma_file.write(">;\n")
        
        magma_file.write('SetNthreads(8);\n')
        
        # TODO add dreg computation here
        
        magma_file.write('SetVerbose("Faugere",2);\n')
        magma_file.write('SetVerbose("FGLM",3);\n')
        
        # GB computation wrt ordering (degrevlex)
        magma_file.write('printf "---------------\\nGB computation\\n---------------\\n";\n')
        magma_file.write('time GB, dregs := GroebnerBasis(I : Al := "Direct", Faugere := true);\n')
        magma_file.write('dregs;\n')
        #magma_file.write('GB;\n')
        magma_file.write('#GB;\n')  # Number of polynomials in GB basis
        magma_file.write('[Degree(f) : f in GB];\n')  # Degrees of polynomials in GB basis
        
        # FGLM basis conversion
        magma_file.write('printf "---------------\\nFGLM computation\\n---------------\\n";\n')
        magma_file.write('I := ChangeOrder(I,"lex");\n')
        magma_file.write('time GB_lex := GroebnerBasis(I: Al := "FGLM");\n')
        #magma_file.write('GB_lex;\n')
        magma_file.write('#GB_lex;\n')  # Number of polynomials in GB_lex basis
        magma_file.write('[Degree(f) : f in GB_lex];\n')  # Degrees of polynomials in GB_lex basis
        #magma_file.write('Degree(gb_lex[#gb_lex]);\n')  # Degree of univariate polynomial
        
        # Variety
        magma_file.write('printf "---------------\\nVARIETY computation\\n---------------\\n";\n')
        magma_file.write('time V := Variety(I);\n')
        magma_file.write('V;\n')
        magma_file.write('#V;\n')
        
        magma_file.write('// exit;\n')     


# ==============================================================================
# Target specific part
# ==============================================================================
def gen_magma_file_tip(tip, eqs, fixed_elems, attacktype, mdir="magma", replace=True):
    # Replace or add counter to filename
    import os
    filename = f"{attacktype}_m{tip.m}r{tip.r}c{tip.c}d{tip.d}_alpha{tip.alpha}_p{tip.p:02x}"
    if replace or not os.path.isfile(mdir + "/" + filename + ".magma"):
        filename = filename + ".magma"
    else:
        counter = 1
        filename = filename + "_{}.magma"
        while os.path.isfile(mdir + "/" + filename.format(counter)):
            counter += 1
        filename = filename.format(counter)
        
    gen_magma_file(filename, eqs, fixed_elems, mdir)


def gen_magma_file_mono(mono, eqs, fixed_elems, attacktype, mdir="magma", replace=True):
    # Replace or add counter to filename
    import os
    filename = f"{attacktype}_t{mono.t}r{mono.r}c{mono.c}d{mono.d}"
    if replace or not os.path.isfile(mdir + "/" + filename + ".magma"):
        filename = filename + ".magma"
    else:
        counter = 1
        filename = filename + "_{}.magma"
        while os.path.isfile(mdir + "/" + filename.format(counter)):
            counter += 1
        filename = filename.format(counter)
        
    gen_magma_file(filename, eqs, fixed_elems, mdir)


# ____                       __  __    _    ____ __  __    _                  _               _ 
#|  _ \ __ _ _ __ ___  ___  |  \/  |  / \  / ___|  \/  |  / \      ___  _   _| |_ _ __  _   _| |_
#| |_) / _` | '__/ __|/ _ \ | |\/| | / _ \| |  _| |\/| | / _ \    / _ \| | | | __| '_ \| | | | __|
#|  __/ (_| | |  \__ \  __/ | |  | |/ ___ \ |_| | |  | |/ ___ \  | (_) | |_| | |_| |_) | |_| | |_
#|_|   \__,_|_|  |___/\___| |_|  |_/_/   \_\____|_|  |_/_/   \_\ \___/ \__,_|\__| .__/ \__,_|\__|
#                                                                               |_|

# ==============================================================================
# General part
# ==============================================================================
load('utils/complexities.sage')

def get_varnames(s):
    pattern = re.compile(r'<([^>]+)>')
    varnames = pattern.search(s).group(1)
    return varnames.split(',')
    
def strlist_to_strlist(s):
    return ''.join(s).replace('\n','').replace('[','').replace(']','').replace(' ','').split(',')


def find_list_start_end(lines, startidx=None, endidx=None):
    assert(startidx is not None or endidx is not None)
    
    if startidx is not None and endidx is None:
        assert(lines[startidx].startswith("["))
        endidx = startidx
        while not lines[endidx].endswith("]\n"):
            endidx += 1
    elif startidx is None and endidx is not None:
        assert(lines[endidx].endswith("]\n"))
        startidx = endidx
        while not lines[startidx].startswith("["):
            startidx -= 1
    else:
        assert(lines[startidx].startswith("["))
        assert(lines[endidx].endswith("]\n"))
        
    return startidx, endidx

def parse_intlist(lines, startidx=None, endidx=None):
    startidx, endidx = find_list_start_end(lines, startidx, endidx)        
    return [int(x) for x in strlist_to_strlist(lines[startidx:endidx+1])], startidx, endidx

def parse_varlist(lines, F, startidx=None, endidx=None):
    startidx, endidx = find_list_start_end(lines, startidx, endidx)
    s = ''.join(lines[startidx:endidx+1])
    pattern = re.compile(r'<\s*([^>]+)\s*>')
    matches = pattern.findall(s)
    parsed_lists = [list(map(F, match.split(','))) for match in matches]
    return parsed_lists, startidx, endidx

def parse_time(s):
    s = s.strip('\n').strip('[r]')[len("Time: "):]
    return float(s)

def parse_magma_output(filename, dirname="magma/results"):
    
    if 'tip' in filename:
        s1, t1, u, m, r, c, d, alpha = tip_parse_filename(filename)
        params = {'s1': s1, 't1': t1, 'u': u, 'n': m, 'r': r, 'c': c, 'd': d, 'alpha': alpha}
    elif 'mono' in filename:
        t, r, c, d = mono_parse_filename(filename)
        params = {'n': t, 'r': r, 'c': c, 'd': d}
    else:
        return None
    
    results = {'p': -1, 'nv': -1, 'ne': -1, 'mdeg': -1, 'degs': [None], 'gelim': {'nv': -1, 'ne': -1, 'mdeg': -1, 'degs': [None], 
               'c': {'FGLM': {'est': (-1,-1), 'exp': (-1,-1)}, 'GB': {'est': (-1,-1), 'exp': (-1,-1)}}},
               'GB': {'t': -1, 'degs': [None], 'ne': -1, 'dregs': [None], 'dexp': None, 'dreg': None, 'dmac': None,'c': {'est': (-1,-1), 'exp': (-1,-1)}}, # Step 1: GB
               'FGLM': {'t': -1, 'degs': [None], 'ne': -1, 'vsdim': None, 'bez': None,'c': {'est': (-1,-1), 'exp': (-1,-1)}},               # Step 2: FGLM
               'ELIM': {'t': -1, 'unideg': None, 'nsols': None, 'sols': [None], 'c': None},     # Step 3: ELIM / Variety
               'params': params, 't': -1, 'mem': -1
              } 
    
    # Get degrees, prime and variables from Magma file
    with open(f"{dirname}/{filename}.magma") as f: # TODO update with correct dir
        lines = f.readlines()
        
        # Line 0: Fixed values, Da2, Db2 # TODO
        results['fixed'] = re.search(r"\[(.*)\];\n", lines[0]).group(1)
        
        # Line 1: Degrees
        degs, _, _ = parse_intlist([lines[1].replace(';','')[3:]], startidx=0, endidx=0)
        results['degs'] = degs
        results['ne'] = len(degs)
        results['mdeg'] = max(degs) if len(degs) > 0 else None
        
        # Line 2: Prime
        results['p'] = int(lines[2][len("p := "):-2])
        
        # Line 3: Ring
        varnames = get_varnames(lines[3])
        results['nv'] = len(varnames)
    
    # Get results from Magma output
    F = GF(results['p'])
    P = PolynomialRing(F, varnames, order='lex')
    gblex = []
    
    with open(f"{dirname}/results/{filename}.txt") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Dreg information
            if line.startswith("GB computation"):
                j = i
                while j >= 0 and not lines[j].startswith("Time"):
                    j -= 1
                
                if j > 0:
                    results['GB']['dreg'] = int(lines[j+1])
                
            # Get GB informations (right before FGLM started)
            if line.startswith("FGLM computation"):
                
                j = i
                while not lines[j].startswith("Time"):
                    j -= 1
                
                # Time, dregs, GB?, neqs, degs
                results['GB']['t'] = parse_time(lines[j])
                
                dregs, _, endidx = parse_intlist(lines, startidx=j+1)
                results['GB']['dregs'] = dregs
                results['GB']['dexp'] = max(dregs) if len(dregs) > 0 else None
                
                results['GB']['degs'], startidx, _ = parse_intlist(lines, endidx=i-2)
                results['GB']['ne'] = int(lines[startidx-1])
                
            if line.startswith("VARIETY computation"):
                j = i
                while not lines[j].startswith("Time"):
                    j -= 1
                
                # GB lex: Time, GBlex?, neqs, degs
                results['FGLM']['t'] = parse_time(lines[j])
                
                results['FGLM']['degs'], startidx, _ = parse_intlist(lines, endidx=i-2)
                results['FGLM']['ne'] = int(lines[startidx-1])
                
                # Check if GBlex was printed
                if lines[j+1].startswith('['):
                    gblex = strlist_to_strlist(lines[j+1:startidx-1])
                    gblex = [P(f) for f in gblex]
                    results['ELIM']['unideg'] = gblex[-1].degree() if gblex[-1].is_univariate() else None
                    
                # Variety: Time, V, nsols
                vcomp = False
                for j in range(i, len(lines)):
                    if lines[j].startswith("Time"):
                        vcomp=True
                        break
                 
                if not vcomp:
                    continue # No variety was computed
                
                results['ELIM']['t'] = parse_time(lines[j])
                results['ELIM']['sols'], _, endidx = parse_varlist(lines, F, startidx=j+1)
                results['ELIM']['nsols'] = int(lines[endidx+1])
            
            if line.startswith("Get rep mat (dim: "):
                results['FGLM']['vsdim'] = int(line[len("Get rep mat (dim: "):-2])
        
            if line.startswith("Quotient dimension: "):
                results['FGLM']['vsdim'] = int(line[len("Quotient dimension: "):-1])
        
        if lines[-3].startswith("Total time: "):
            #print(lines[-3])
            # Total time: 341.320 [340.693] seconds, Total memory usage: 2666.97MB
            numbers = re.findall(r'\d+\.\d+|\d+', lines[-3])
            numbers = [float(num) if '.' in num else int(num) for num in numbers]
            #print(numbers)
            results['t'] = numbers[0]
            results['mem'] = numbers[-1]
        
    # Number of variables / equations after performing Gaussian elimination
    num_lineqs = sum([d == 1 for d in results['degs']])
    nv = results['nv'] - num_lineqs
    ne = results['ne'] - num_lineqs
    degs = [d for d in results['degs'] if d != 1]
    results['gelim']['nv'] = nv
    results['gelim']['ne'] = ne
    results['gelim']['degs'] = degs
    results['gelim']['mdeg'] = max(degs)
    
    # Not affected by removing linear polynomials
    macaulay = lambda degs: sum(map(lambda x: x-1, degs)) + 1
    bezout = lambda degs: prod(degs)

    # ----- GB complexity -----
    dmac = macaulay(degs)
    results['GB']['dmac'] = dmac
    
    # Estimated GB complexity (with a priori Gaussian elimination)
    k, d = gb_compl_kd(dmac, nv, ne)
    results['gelim']['c']['GB']['est'] = (k,d)
    
    # Estimated GB complexity
    k, d = gb_compl_kd(dmac, results['nv'], results['ne'])
    results['GB']['c']['est'] = (k,d)
    
    # Experimental GB complexity (with a priori Gaussian elimination, assuming dreg is the same)
    k, d = gb_compl_kd(results['GB']['dreg'], nv, ne)
    results['gelim']['c']['GB']['exp'] = (k,d)
    
    # Experimental GB complexity
    k, d = gb_compl_kd(results['GB']['dreg'], results['nv'], results['ne'])
    results['GB']['c']['exp'] = (k,d)
    
    # ----- FGLM complexity -----
    bez = bezout(degs)
    results['FGLM']['bez'] = bez
    
    # Estimated FGLM complexity (with a priori Gaussian elimination)              
    k, d = fglm_compl_kd(bez, nv, ne)
    results['gelim']['c']['FGLM']['est'] = (k,d)
    
    # Estimated FGLM complexity
    k, d = fglm_compl_kd(bez, results['nv'], results['ne'])
    results['FGLM']['c']['est'] = (k,d)
    
    # Experimental FGLM complexity (with a priori Gaussian elimination)    
    k, d = fglm_compl_kd(results['FGLM']['vsdim'], nv, ne)
    results['gelim']['c']['FGLM']['exp'] = (k,d)
    
    # Experimental FGLM complexity
    k, d = fglm_compl_kd(results['FGLM']['vsdim'], results['nv'], results['ne'])
    results['FGLM']['c']['exp'] = (k,d)
             
    return results, gblex, P


# ==============================================================================
# Target specific part
# ==============================================================================
def tip_parse_filename(filename):
    # Regular expression with named groups
    pattern = re.compile(r's1_(?P<s1>\d+)_t1_(?P<t1>\d+)(_u_(?P<u>\d+))?_m(?P<m>\d+)r(?P<r>\d+)c(?P<c>\d+)d(?P<d>\d+)_alpha(?P<alpha>\d+)')

    # Search for matches
    match = pattern.search(filename)

    # Extract values using named groups
    s1 = int(match.group('s1')) if match else None
    t1 = int(match.group('t1')) if match else None
    u = int(match.group('u')) if match and match.group('u') else None
    m = int(match.group('m')) if match else None
    r = int(match.group('r')) if match else None
    c = int(match.group('c')) if match else None
    d = int(match.group('d')) if match else None
    alpha = int(match.group('alpha')) if match else None
    
    return s1, t1, u, m, r, c, d, alpha


def mono_parse_filename(filename):
    # Regular expression with named groups
    pattern = re.compile(r'_t(?P<t>\d+)r(?P<r>\d+)c(?P<c>\d+)d(?P<d>\d+)')

    # Search for matches
    match = pattern.search(filename)

    # Extract values using named groups
    t = int(match.group('t')) if match else None
    r = int(match.group('r')) if match else None
    c = int(match.group('c')) if match else None
    d = int(match.group('d')) if match else None
    
    return t, r, c, d