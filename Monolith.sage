
# The defining, first column of the (circulant) MDS matrix.
# See, for example:
# https://extgit.iaik.tugraz.at/krypto/zkfriendlyhashzoo/-/blob/master/plain_impls/src/monolith_31/mds_16_generic.rs?ref_type=heads
# https://extgit.iaik.tugraz.at/krypto/zkfriendlyhashzoo/-/blob/master/plain_impls/src/monolith_31/mds_24_generic.rs?ref_type=heads
_MDS8 =  [23, 8, 13, 10, 7, 6, 21, 8]
_MDS12 = [7, 23, 8, 26, 13, 10, 9, 7, 6, 22, 21, 8]
_MDS16 = [61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689, 26798, 17845]
_MDS32 = [0x536C316, 0x1DD20A84, 0x43E26541, 0x52B22B8D, 0x37DABDF0, 0x540EC006, 0x3015718D, 0x5A99E14C, 
          0x23637285, 0x4C8A2F76, 0x5DEC4E6E, 0x374EE8D6, 0x27EDA4D8, 0x665D30D3, 0x32E44597, 0x43C7E2B3, 
          0x67C4C603, 0x78A8631F, 0x452F77E3, 0x39F03DF, 0x743DBFE0, 0x4DA05A48, 0x5F027940, 0x8293632,
          0x50F2C76A, 0x7B773729, 0x577DE8B0, 0x73B1EAC6, 0x58DA7D29, 0x67AA4375,0xDBA9E33, 0x2655E5A1]

class Monolith:
    
    # -----------------------------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------------------------
    
    def __init__(self, n=64, eta=32, t=12, u=None, sis=[8, 8, 8, 8, 8, 8, 8, 8], N=1, r=None, c=None, d=None, debug=False):
        self.debug = debug
        # construct prime field
        if eta is None:
            # Mersenne prime
            p = 2^n - 1
        else:
            # Goldilocks prime
            p = 2^n - 2^eta + 1
        self.name = f"Mono{n}"
        assert(is_prime(p))
        self.n = n
        self.eta = eta
        self.p = p
        self.F = GF(p)
        # inputs and outputs to permutation in F^t
        # 8 for compression function, 12 for sponge usage
        # assert(t in [8, 12])
        self.t = t
        # Decomposition information for BARS layer, s is the size of the splits
        # Bars(x0,...,xt-1) = Bar(x0) || Bar(x1) || .. || Bar(xu-1) || xu || ... || xt-1
        if u is None:
            self.u = t
        else:
            assert(u <= t)
            self.u = u
        self.sis = sis
        self.m = len(sis)
        # Number of rounds
        self.N = N
        
        # Rate, capacity and digest length, if used in Sponge hashing mode
        self.r = r
        self.c = c
        self.d = d
        self.c_value = 0 # TODO check with safe mode
        
        # BIT MANIPULATION FUNCTIONS
        self.bitset = lambda x, b : ((x >> (b - 1)) > 0)
        self.bitinv = lambda x, s : ((x ^^ (2**s - 1)) & (2**s - 1))
        self.rotl = lambda x, s, b: ((x << b) & (2**s - 1)) | (x >> (s - b))
        self.dmsb = lambda x, s, b: (x & ((2**s - 1) >> b)) # Discard b MSBs for x (which is of bitsize s)
        
        # CONVERSION HELPER FUNCTIONS: FINITE FIELD ELEMENT <-> INTEGER
        self.i2f = lambda x : self.F(x)
        self.f2i = lambda x : Integer(x)
            
        # STRING REPRESENTATION FOR DECOMPOSED VALUES
        self.l2s = lambda xis : [f"{xi:0{si}b}" for xi,si in zip(xis,self.sis)]
        
        # MDS matrix and round constants
        self.M = self.get_m_concrete()
        self.M_inv = self.M.inverse()
        self.rcons = self.get_rcons()
    
    
    def __repr__(self):
        repr_string = f"{self.name} permutation\n"
        repr_string += f"  p = 2^{self.n} - 2^{self.eta} + 1 = {self.p}\n"
        repr_string += f"  t = {self.t}, u = {self.u}, m = {self.m}, sis = {self.sis}\n"
        repr_string += f"  Number of rounds: {self.N}\n"
        return repr_string
    
    
    # -----------------------------------------------------------------------------------------
    # BARS
    # -----------------------------------------------------------------------------------------
    
    def decompose(self, x):
        assert(x in ZZ)
        
        x_ = x
        xis = []
        for si in self.sis[::-1]:
            xis.append(x_ & (2**si - 1))
            x_ = x_ >> si # next chunk
        xis = xis[::-1] # chunks[0] contains most significant bits
        
        if self.debug:
            print(f"DECOMPOSE: {x} -> {self.l2s(xis)}")
        
        assert(all([xi < 2**si for xi,si in zip(xis, self.sis)]))
        assert(int(''.join(self.l2s(xis)),2) == x)
        return xis
    
    def compose(self, xis):
        assert(all([xi in ZZ for xi in xis]))
        assert(all([xi < 2**si for xi,si in zip(xis, self.sis)]))
        
        x = xis[0]
        for xi,si in zip(xis[1:], self.sis[1:]):
            x = x << si  
            x += xi
        
        if self.debug:
            print(f"COMPOSE: {self.l2s(xis)} -> {x}")
        
        assert(x < self.p)
        assert(int(''.join(self.l2s(xis)),2) == x)
        return x
    
    def sbox(self, x, s):
        # Internal S-Box (Chi-function): S: F_2^s -> F_2^s 
        assert(x in ZZ)
        assert(x < 2**s) # ensure that x in GF(2^s)
        
        t1 = self.rotl(self.bitinv(x,s), s, 1)
        t2 = self.rotl(x, s, 2)
        
        if gcd(s,2) == 1 or s == 2: # s=7, and workaround for tiny instances with s=2
            z = x ^^ (t1 & t2)
        elif gcd(s,3) == 1: # s = 8
            t3 = self.rotl(x, s, 3)
            z = x ^^ (t1 & t2 & t3)
        else:
            raise Exception(f"S: F_2^{s} -> F_2^{s} not implemented.")
        
        z = self.rotl(z, s, 1)
        assert(z < 2**s) # ensure that z in GF(2^s)
        
        if self.debug:
            print(f"S: F_2^{s} -> F_2^{s}. S({x:0{s}b}) = {z:0{s}b}")
        
        return z
    
    
    # Apply BAR to element of GF(p): BAR(x) = COMPOSE(S(DECOMPOSE(x)))
    def bar(self, x, inv=False):
        assert(x in self.F)
        assert(inv==False)
        
        xs = self.decompose(self.f2i(x))
        ys = []
        
        for xi, si in zip(xs, self.sis):
            ys.append(self.sbox(xi, si))
        y = self.i2f(self.compose(ys))
        
        assert(y in self.F)
        return y
    
    
    def bars(self, state, inv=False): # TODO implement inversion?
        state_cpy = copy(state)
        for i in range(self.u):
            state_cpy[i] = self.bar(state[i], inv)
        return state_cpy
    
    
    # -----------------------------------------------------------------------------------------
    # BRICKS
    # -----------------------------------------------------------------------------------------
    
    def bricks(self, state, inv=False):
        out = []
        out.append(state[0])
        for i in range(1, self.t):
            out.append((state[i] - out[i-1]^2) if inv else (state[i] + state[i-1]^2))
        assert(len(out) == self.t)

        return vector(out)
    
    # -----------------------------------------------------------------------------------------
    # CONCRETE
    # -----------------------------------------------------------------------------------------
    
    # Get matrix for CONCRETE
    def get_m_concrete(self):
        from scipy.linalg import circulant
        
        if self.p == 2**64 - 2**32 + 1:
            # Double-checked MDS matrices with reference implementation
            # https://extgit.iaik.tugraz.at/krypto/zkfriendlyhashzoo/-/tree/master/plain_impls/src/monolith_64?ref_type=heads
            if self.t == 8: # 2:1 compression Monolith-64
                return Matrix(self.F, circulant(_MDS8)).transpose()
            elif self.t == 12: # Sponge Monolith-64
                return Matrix(self.F, circulant(_MDS12)).transpose()
            else:
                assert(False)
        elif self.p == 2**31 - 1:
            # Double-checked MDS matrices with reference implementation
            # https://extgit.iaik.tugraz.at/krypto/zkfriendlyhashzoo/-/tree/master/plain_impls/src/monolith_31?ref_type=heads
            if self.t == 16: # 2:1 compression Monolith-31
                return Matrix(self.F, circulant(_MDS16))
            elif self.t == 24: # Sponge Monolith-31
                M = Matrix(self.F, circulant(_MDS32)).transpose()
                return M.submatrix(row=0,col=0,nrows=self.t,ncols=self.t)
            else:
                assert(False)
        else:
            assert(False)
            #return Matrix(self.F, circulant([2] + [1]*(self.t-1)))
            
    
    # CONCRETE layer with round constant addition
    def concrete(self, state, inv=False):
        return self.M_inv * state if inv else self.M * state
    
    # -----------------------------------------------------------------------------------------
    # ROUND CONSTANTS
    # -----------------------------------------------------------------------------------------
    
    # Get round constants 
    def get_rcons(self):        
        from Crypto.Hash import SHAKE128
        shake = SHAKE128.new()
        
        ### SHAKE128 INIT with Monolith<t><N><p in little endian><partition>
        num_bytes_p = ceil(ceil(log(self.p, 2)) / 8.0)
        init_data = b'Monolith' + int(self.t).to_bytes(1, 'little') + int(self.N).to_bytes(1, 'little') 
        init_data += int(self.p).to_bytes(num_bytes_p, 'little') 
        init_data += b''.join([int(si).to_bytes(1, 'little') for si in self.sis[::-1]])
        shake.update(init_data)
        
        def field_element_from_shake(shake):
            bitlen = ceil(log(self.p, 2))
            byte = ceil(bitlen / 8)
            word = ceil(byte / 8)

            while True:
                buf = shake.read(int(byte))
                cons = int.from_bytes(buf, "little")
                if cons < self.p:
                    return self.i2f(cons)
        
        # rcons zero in last round
        # t*(N-1) round constants, t*0 at end            
        rcons = []
        while len(rcons) != self.t * (self.N - 1):
            rcons.append(field_element_from_shake(shake))
                
        # split flat list of round constants into chunks of size t
        rcons = [rcons[i:i + self.t] for i in range(0, len(rcons), self.t)]
        rcons += [[self.F.zero()] * self.t]
        
        assert(all([rcon <= self.p - 1 for rcon in flatten(rcons)])) # Comparison in GF(p), where p == 0
        for r in range(len(rcons)):
            rcons[r] = vector(self.F, rcons[r])
        
        return rcons
    
    
    def add_rcons(self, state, r, inv=False):
        return state - self.rcons[r] if inv else state + self.rcons[r]
    
    
    # -----------------------------------------------------------------------------------------
    # PERMUTATION
    # -----------------------------------------------------------------------------------------
        
    def linear_layer(self, state, r, inv=False):
        if inv:   
            state = self.add_rcons(state, r, inv)
            state = self.concrete(state, inv)
        else:
            state = self.concrete(state, inv)
            state = self.add_rcons(state, r, inv)
        return state
    
    def nonlinear_layer(self, state, inv=False):
        if inv:
            state = self.bricks(state, inv)
            state = self.bars(state, inv)
        else:
            state = self.bars(state, inv)
            state = self.bricks(state, inv)
        return state
    
    def __call__(self, message):
        assert(len(message) == self.r) # Only rate part can be chosen, capacity is fixed
        input_state = vector(self.F, message + [self.c_value]*self.c) # Append capacity
        result = self.eval_with_intermediate_states(input_state)
        hash_value = list(result[-1])[:self.d]
        return hash_value
        
    def eval_with_intermediate_states(self, input_state):
        # input_state is a list of input values 
        if len(input_state) == self.r:
            input_state = input_state + [self.c_value]*self.c
        input_state = vector(self.F, input_state)
        assert(len(input_state) == self.t)
        
        result = [input_state]

        for r in range(self.N):
            state = result[-1]
            
            # Initial application of CONCRETE layer
            if r == 0:
                state = self.concrete(state)
                
            state = self.nonlinear_layer(state)
            state = self.linear_layer(state, r)
            
            result.append(state)
        return result

    
# -----------------------------------------------------------------------------------------
# Instances
# -----------------------------------------------------------------------------------------

Mono64_sponge = Monolith(n=64, eta=32, t=12, u=4, sis=[8, 8, 8, 8, 8, 8, 8, 8], r=8, c=4, d=4, N=6)
Mono31_sponge = Monolith(n=31, eta=None, t=24, u=8, sis=[7, 8, 8, 8], r=16, c=8, d=8, N=6)

Mono64_compression = Monolith(n=64, eta=32, t=8, u=4, sis=[8, 8, 8, 8, 8, 8, 8, 8], N=6)
Mono31_compression = Monolith(n=31, eta=None, t=16, u=8, sis=[7, 8, 8, 8], N=6)


TV_monolith = {
    '64': {
        'sponge' : [
            {
                'in':  list(IntegerRange(0,Mono64_sponge.t)),
                'out' : [5867581605548782913, 588867029099903233, 6043817495575026667, 805786589926590032, 9919982299747097782, 6718641691835914685, 7951881005429661950, 15453177927755089358, 974633365445157727, 9654662171963364206, 6281307445101925412, 13745376999934453119]
            }
        ],
        'compression' : [
            {
                'in' : list(IntegerRange(0,Mono64_compression.t)),
                'out' : [3656442354255169651, 1088199316401146975, 22941152274975507, 14434181924633355796, 6981961052218049719, 16492720827407246378, 17986182688944525029, 9161400698613172623]
            }
        ]
    },
    '31': {
        'sponge' : [
            {
                'in' : list(IntegerRange(0,Mono31_sponge.t)),
                'out' : [2067773075, 1832201932, 1944824478, 1823377759, 1441396277, 2131077448, 2132180368, 1432941899, 1347592327, 1652902071, 1809291778, 1684517779, 785982444, 1037200378, 1316286130, 1391154514, 1760346031, 1412575993, 2108791223, 1657735769, 219740691, 1165267731, 505815021, 2080295871]
            }
        ],
        'compression' : [
            {
                'in' : list(IntegerRange(0,Mono31_compression.t)),
                'out' : [609156607, 290107110, 1900746598, 1734707571, 2050994835, 1648553244, 1307647296, 1941164548, 1707113065, 1477714255, 1170160793, 93800695, 769879348, 375548503, 1989726444, 1349325635]
            }
        ]
    }
}


# Test PERMUTATION
for test_vector in TV_monolith['64']['sponge']:
    perm_in = test_vector['in']
    perm_out = Mono64_sponge.eval_with_intermediate_states(perm_in)[-1]
    assert(list(perm_out) == test_vector['out'])

for test_vector in TV_monolith['64']['compression']:
    perm_in = test_vector['in']
    perm_out = Mono64_compression.eval_with_intermediate_states(perm_in)[-1]
    assert(list(perm_out) == test_vector['out'])
    
for test_vector in TV_monolith['31']['compression']:
    perm_in = test_vector['in']
    perm_out = Mono31_compression.eval_with_intermediate_states(perm_in)[-1]
    assert(list(perm_out) == test_vector['out'])
    
for test_vector in TV_monolith['31']['sponge']:
    perm_in = test_vector['in']
    perm_out = Mono31_sponge.eval_with_intermediate_states(perm_in)[-1]
    assert(list(perm_out) == test_vector['out'])