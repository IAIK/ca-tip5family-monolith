# --------------------------------------------------------------------------------------------
# Constants from Tip5 Rust implementation
# https://github.com/Neptune-Crypto/twenty-first/blob/master/twenty-first/src/math/tip5.rs
# --------------------------------------------------------------------------------------------

from hashlib import sha256
from blake3 import blake3

SHA256 = lambda x: sha256(x).digest()
Blake3 = lambda x, n: blake3(x).digest(n)

LOOKUP_TABLE = [
    0, 7, 26, 63, 124, 215, 85, 254, 214, 228, 45, 185, 140, 173, 33, 240, 29, 177, 176, 32, 8,
    110, 87, 202, 204, 99, 150, 106, 230, 14, 235, 128, 213, 239, 212, 138, 23, 130, 208, 6, 44,
    71, 93, 116, 146, 189, 251, 81, 199, 97, 38, 28, 73, 179, 95, 84, 152, 48, 35, 119, 49, 88,
    242, 3, 148, 169, 72, 120, 62, 161, 166, 83, 175, 191, 137, 19, 100, 129, 112, 55, 221, 102,
    218, 61, 151, 237, 68, 164, 17, 147, 46, 234, 203, 216, 22, 141, 65, 57, 123, 12, 244, 54, 219,
    231, 96, 77, 180, 154, 5, 253, 133, 165, 98, 195, 205, 134, 245, 30, 9, 188, 59, 142, 186, 197,
    181, 144, 92, 31, 224, 163, 111, 74, 58, 69, 113, 196, 67, 246, 225, 10, 121, 50, 60, 157, 90,
    122, 2, 250, 101, 75, 178, 159, 24, 36, 201, 11, 243, 132, 198, 190, 114, 233, 39, 52, 21, 209,
    108, 238, 91, 187, 18, 104, 194, 37, 153, 34, 200, 143, 126, 155, 236, 118, 64, 80, 172, 89,
    94, 193, 135, 183, 86, 107, 252, 13, 167, 206, 136, 220, 207, 103, 171, 160, 76, 182, 227, 217,
    158, 56, 174, 4, 66, 109, 139, 162, 184, 211, 249, 47, 125, 232, 117, 43, 16, 42, 127, 20, 241,
    25, 149, 105, 156, 51, 53, 168, 145, 247, 223, 79, 78, 226, 15, 222, 82, 115, 70, 210, 27, 41,
    1, 170, 40, 131, 192, 229, 248, 255,
]


# NUM_ROUNDS * STATESIZE = 5 * 16
ROUND_CONSTANTS = [
    13630775303355457758,
    16896927574093233874,
    10379449653650130495,
    1965408364413093495,
    15232538947090185111,
    15892634398091747074,
    3989134140024871768,
    2851411912127730865,
    8709136439293758776,
    3694858669662939734,
    12692440244315327141,
    10722316166358076749,
    12745429320441639448,
    17932424223723990421,
    7558102534867937463,
    15551047435855531404,
    17532528648579384106,
    5216785850422679555,
    15418071332095031847,
    11921929762955146258,
    9738718993677019874,
    3464580399432997147,
    13408434769117164050,
    264428218649616431,
    4436247869008081381,
    4063129435850804221,
    2865073155741120117,
    5749834437609765994,
    6804196764189408435,
    17060469201292988508,
    9475383556737206708,
    12876344085611465020,
    13835756199368269249,
    1648753455944344172,
    9836124473569258483,
    12867641597107932229,
    11254152636692960595,
    16550832737139861108,
    11861573970480733262,
    1256660473588673495,
    13879506000676455136,
    10564103842682358721,
    16142842524796397521,
    3287098591948630584,
    685911471061284805,
    5285298776918878023,
    18310953571768047354,
    3142266350630002035,
    549990724933663297,
    4901984846118077401,
    11458643033696775769,
    8706785264119212710,
    12521758138015724072,
    11877914062416978196,
    11333318251134523752,
    3933899631278608623,
    16635128972021157924,
    10291337173108950450,
    4142107155024199350,
    16973934533787743537,
    11068111539125175221,
    17546769694830203606,
    5315217744825068993,
    4609594252909613081,
    3350107164315270407,
    17715942834299349177,
    9600609149219873996,
    12894357635820003949,
    4597649658040514631,
    7735563950920491847,
    1663379455870887181,
    13889298103638829706,
    7375530351220884434,
    3502022433285269151,
    9231805330431056952,
    9252272755288523725,
    10014268662326746219,
    15565031632950843234,
    1209725273521819323,
    6024642864597845108,
]

def Gen_RoundConstant(p, m=16, N=5):
    '''
    checked with constants in [1]
    '''
    RC = [[0 for _ in range(m)] for _ in range(N)]
    Fp = FiniteField(p)

    mul = Fp(2**64)
    mul_inverse = mul^-1

    for r in range(N):
        for i in range(m):
            # get byte str first
            i_value = int(i + r * m)
            bytes_val = i_value.to_bytes(1, 'big')
            # print(bytes_val)
            seed_str = "{}".format("Tip5")
            # concat then hash
            byte_str = Blake3(bytes(seed_str, "ascii") + bytes_val, m)

            # process bytes into integer
            integer = Fp(sum((2**8)**j * byte_str[j] for j in range(m)))
            integer *= mul_inverse
            RC[r][i] = integer
        RC[r] = vector(Fp, RC[r])

    return RC


MDS_MATRIX_FIRST_COLUMN = [
    61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689,
    26798, 17845,
]


# --------------------------------------------------------------------------------------------
# SAGE implementation Tip5 (and Tip4)
# --------------------------------------------------------------------------------------------

from sage.crypto.sbox import SBox
from scipy.linalg import circulant
from random import randint

class Tip5:
    
    def __init__(self, p=2**64 - 2**32 + 1, state_size=16, num_rounds=5, alpha=7, name="TIP", num_split_and_lookup=4, rate=10, capacity=6, digest_length=5, initial_capacity_value=1, R=2**64):
        assert(ceil(log(p,2)) == 64)
        self.p = p # 0xffffffff00000001 = 2**64 - 2**32 + 1
        self.p_orig = 0xffffffff00000001
        self.field = GF(self.p) 
        self.to_field   = lambda x : self.field(x)
        self.from_field = lambda x : Integer(x)
        self.random_element = lambda : self.to_field(randint(0,min(self.p_orig,self.p)-1)) # Make sure that random inputs work for S S-Box
        
        max_state_size, max_num_rounds = 16, 5
        assert(state_size <= max_state_size and num_rounds <= max_num_rounds)
        self.m = state_size
        self.N = num_rounds
        self.s = num_split_and_lookup
        
        # Components for nonlinear layer: S
        self.lut = SBox([x for x in LOOKUP_TABLE]) # Lookup table defined over GF(2**8+1), elements interpreted as integers in LUT
        self.lut_inv = self.lut.inverse()
        self.R = self.to_field(R)
        self.R_inv = self.R**(-1)
        
        # Components for nonlinear layer: T
        assert(gcd(alpha,self.p - 1)==1)
        self.alpha = alpha
        self.alpha_inv = inverse_mod(self.alpha, self.p - 1)
        
        # Components for linear layer
        self.mds_matrix = Matrix(self.field, nrows=self.m, ncols=self.m, entries=circulant(MDS_MATRIX_FIRST_COLUMN[:self.m]))
        self.mds_matrix_inv = self.mds_matrix.inverse()
        #self.rcons = [vector(self.field, ROUND_CONSTANTS[r*self.m:(r+1)*self.m]) for r in range(self.N)]
        self.rcons = Gen_RoundConstant(self.p, self.m, self.N)
        
        # Hash mode specifics
        assert(self.m == rate + capacity)
        self.r = rate
        self.c = capacity
        assert(self.m >= digest_length)
        self.d = digest_length
        self.c_value = self.to_field(initial_capacity_value)
        
        self.name = name
    
    def __str__(self):
        result = f"{self.name} over {self.field}\n"
        result += f"Rate = {self.r}, Capacity = {self.c}, Digest length = {self.d}, Number of split and look S-Boxes = {self.s}, Alpha = {self.alpha}\nNumber of rounds = {self.N}, Initial capacity value = {self.c_value}, R for Montgomery = {self.R}\n"
        return result
    
    # --------------------------------------------------------------------------------
    # Call hash function
    # --------------------------------------------------------------------------------    
    def __call__(self, message):
        assert(len(message) == self.r) # Only rate part can be chosen, capacity is fixed
        input_state = vector(self.field, message + [self.c_value]*self.c) # Append capacity
        result = self.eval_with_intermediate_states(input_state)
        hash_value = list(result[-1])[:self.d]
        return hash_value
    
    def eval_with_intermediate_states(self, input_state, N=None):
        if len(input_state) == self.r:
            input_state = input_state + [self.c_value]*self.c
        input_state = vector(self.field, input_state)
        assert(len(input_state) == self.m)
        # If no custom round number N is given, use predefined one
        if N is None:
            N = self.N
        
        result = [input_state]
        for r in range(N):
            state = result[-1]
            state = self.nonlinear_layer(state)
            state = self.linear_layer(state, r)
            result.append(state)
        return result
    
    def eval_inv_with_intermediate_states(self, output_state):
        assert(len(output_state) == self.m) # Inversion only possible if output rate and capacity are known
        result = [vector(self.field, output_state)]
        for r in range(self.N-1, -1, -1):
            state = result[-1]
            state = self.linear_layer(state, r, inv=True)
            state = self.nonlinear_layer(state, inv=True)
            result.append(state)
        return result
    
    # --------------------------------------------------------------------------------
    # Non-linear layer: Substitution S (with composition and decomposition) and T (power map)
    # --------------------------------------------------------------------------------
    def decompose(self, x):
        x = self.from_field(x)
        parts = []
        for si in [8]*8: # 64bit = 8x8 bits
            parts.append(x & (2**si - 1))
            x = x >> si # next chunk
        return parts[::-1] # parts[0] contains most significant 8 bits

    def compose(self, parts):
        x = parts[0]
        for xi,si in zip(parts[1:], [8]*7): # 64bit = 8bits + 7x8 bits
            x = x << si  
            x += xi
        assert(x < self.p)
        return self.to_field(x)

    def S(self, x, inv=False):
        x = self.R * x  # Canonical Representation -> Montgomery Form
        parts = self.decompose(x)
        sub_parts = [self.lut_inv(xi) if inv else self.lut(xi) for xi in parts]
        y = self.compose(sub_parts)
        y = self.R_inv * y  # Montgomery Form -> Canonical Representation
        return y
    
    def T(self, x, inv=False):
        return x**self.alpha_inv if inv else x**self.alpha
    
    def nonlinear_layer(self, state, inv=False):
        sub_state = copy(state)
        for i in range(self.m):
            if i < self.s:
                sub_state[i] = self.S(sub_state[i], inv)
            else:
                sub_state[i] = self.T(sub_state[i], inv)
        return sub_state
    
    # --------------------------------------------------------------------------------
    # Linear layer: MDS matrix multiplication and round constant addition
    # --------------------------------------------------------------------------------
    def multiply_mds_matrix(self, state, inv=False):
        return self.mds_matrix_inv * state if inv else self.mds_matrix * state
    
    def add_rcons(self, state, r, inv=False):
        #return state - self.rcons.column(r) if inv else state + self.rcons.column(r)
        return state - self.rcons[r] if inv else state + self.rcons[r]
    
    def linear_layer(self, state, r, inv=False):
        if inv:   
            state = self.add_rcons(state, r, inv)
            state = self.multiply_mds_matrix(state, inv)
        else:
            state = self.multiply_mds_matrix(state, inv)
            state = self.add_rcons(state, r, inv)
        return state


# --------------------------------------------------------------------------------------------
# Test vectors created by Rust implementation: hash10_test_vectors
# https://github.com/Neptune-Crypto/twenty-first/blob/master/twenty-first/src/math/tip5.rs
# --------------------------------------------------------------------------------------------

TV_Tip5 = [
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [941080798860502477, 5295886365985465639, 14728839126885177993, 10358449902914633406, 14220746792122877272]],
    [[941080798860502477, 5295886365985465639, 14728839126885177993, 10358449902914633406, 14220746792122877272, 0, 0, 0, 0, 0], [15888421881075650037, 8699648354187865464, 6719068786850902915, 16188941274693647820, 4768361305800190493]],
    [[941080798860502477, 15888421881075650037, 8699648354187865464, 6719068786850902915, 16188941274693647820, 4768361305800190493, 0, 0, 0, 0], [11494362724359741120, 2984169814429715553, 11021746812971026026, 5102281498552384717, 5023112854146751042]],
    [[941080798860502477, 15888421881075650037, 11494362724359741120, 2984169814429715553, 11021746812971026026, 5102281498552384717, 5023112854146751042, 0, 0, 0], [627201255727529993, 2530132417472465719, 15134374672529870482, 10586143339158028166, 13810271029904013559]],
    [[941080798860502477, 15888421881075650037, 11494362724359741120, 627201255727529993, 2530132417472465719, 15134374672529870482, 10586143339158028166, 13810271029904013559, 0, 0], [4790238723037855394, 13717377209729127271, 8994982932799814404, 18004412270774820131, 5877166878145340765]],
    [[941080798860502477, 15888421881075650037, 11494362724359741120, 627201255727529993, 4790238723037855394, 13717377209729127271, 8994982932799814404, 18004412270774820131, 5877166878145340765, 0], [16959020643814878453, 12118009629857908438, 10239930869937551135, 6889489196156760098, 5774309862903741805]],
    [[941080798860502477, 15888421881075650037, 11494362724359741120, 627201255727529993, 4790238723037855394, 16959020643814878453, 12118009629857908438, 10239930869937551135, 6889489196156760098, 5774309862903741805], [10869784347448351760, 1853783032222938415, 6856460589287344822, 17178399545409290325, 7650660984651717733]]
]


# --------------------------------------------------------------------------------
# Tip5
# --------------------------------------------------------------------------------
tip5 = Tip5(name="Tip5")
#print("tip5 = ")
#print(tip5)

# Check test vectors for Tip5
for tv in TV_Tip5:
    assert(tip5(tv[0]) == tv[1])
    
results = tip5.eval_with_intermediate_states(TV_Tip5[0][0] + [tip5.c_value]*tip5.c)
assert(tip5.eval_inv_with_intermediate_states(results[-1])[-1] == results[0])

# --------------------------------------------------------------------------------
# Tip4
# --------------------------------------------------------------------------------
tip4 = Tip5(name="TIP4", rate=12, capacity=4, digest_length=4)
#print("tip4 = ")
#print(tip4)

# --------------------------------------------------------------------------------
# Tip4'
# --------------------------------------------------------------------------------
tip4p = Tip5(name="TIP4p", state_size=12, rate=8, capacity=4, digest_length=4)
#print("tip4p = ")
#print(tip4p)

# --------------------------------------------------------------------------------
# TipToy
# --------------------------------------------------------------------------------
tiptoy = Tip5(name="TIPTOY", p=2**64-0x3b, num_rounds=3, state_size=6, rate=4, capacity=2, digest_length=2, alpha=3, num_split_and_lookup=2, R=1, initial_capacity_value=0)
#print("tiptoy = ")
#print(tiptoy)

# --------------------------------------------------------------------------------
# MiniTipToy (solveable in sage)
# --------------------------------------------------------------------------------
minitiptoy = Tip5(name="TIPTOYMINI", p=2**64-0x3b, num_rounds=3, state_size=3, rate=2, capacity=1, digest_length=1, alpha=3, num_split_and_lookup=1, R=1, initial_capacity_value=0)
#print("minitiptoy = ")
#print(minitiptoy)

# --------------------------------------------------------------------------------
# Tip2 (Tip5 but smaller alpha, halved state size and adapted rate + capacity)
# --------------------------------------------------------------------------------
tip2 = Tip5(name="tip2", p=2**64-0x3b, num_rounds=3, state_size=8, rate=5, capacity=3, digest_length=2, alpha=3, num_split_and_lookup=2, R=1, initial_capacity_value=0)