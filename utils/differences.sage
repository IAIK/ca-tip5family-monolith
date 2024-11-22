# Helper functions
from itertools import product
from random import choices, randint


# DDT for Lookuptable "U" (8-bit) used in S-box S, with the additional restriction that U(x_1) < U(x_2) for x_1 < x_2.
# This is a heavily simplified version of the presented algorithm, but enough to show that attacks are practical.
TIP5_LOOKUP_TABLE = [
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

DDT_Uspec = [[0 for _ in range(2**8)] for _ in range(2**8)]
for dx in range(1,2**8):
    for x in range(2**8-dx):
        dy = TIP5_LOOKUP_TABLE[x + dx] - TIP5_LOOKUP_TABLE[x]
        if dy > 0: # U(x) < U(x + dx) for x < x + dx 
            DDT_Uspec[dx][dy] += 1
DDT_Uspec[0][0] = 2**8

def din_to_dout_S(tip, din_S, u=8, debug=False):
    """
    For a given input difference `din_S` to an S-box, computes the valid output differences `dout_S`.
    Note that multiplication with tip.R (Montgomery constant) is considered outside of this function.
    
    Additional constraints: 
    - For a tuple (xi,xi+di) with xi < xi+di, it shall hold that U(xi) < U(xi + di).
    - Only activate u least most significant parts.
    
    Args:
    - tip (object): A member of the Tip5 family, used to decompose and compose differences.
    - din_S (int): The input difference to the S-box, din_S ~ [d1,d2,...,d8].
    - u (int, optional): The number of least significant bits to be considered active (between 0 and 8). Defaults to 8 (all bits active).
    - debug (bool, optional): If set to True, prints intermediate steps for debugging purposes. Defaults to False.
    
    Returns:
    - (din_S, dout_S) tuple:
        - din_S (int or None): The input difference after applying constraints (only the least significant `u` bits are active), or `None` if no valid IO difference pair can be found.
        - dout_S (int or None): One randomly chosen output difference for the given input difference, or `None` if no valid IO difference pair can be found.
    """
    assert(u >= 0 and u <= 8)
    
    # Only non-zero difference ok here
    if din_S == 0:
        return None, None
    
    dins_U = tip.decompose(din_S) # [d1,d2,...,d8]
    if debug:
        print(f"{din_S} = {dins_U}")
    # only activate last u parts
    for i in range(8-u):
        dins_U[i] = 0
    din_S = tip.compose(dins_U) 
    
    douts_U = []
    for din in dins_U:
        # Check if for given input difference, there is an output difference such that
        # U(x) < U(x + din) for x < x + din 
        if sum(DDT_Uspec[din]) == 0:
            return None, None
        else:
            #print([i for i in range(len(DDT_Uspec[din])) if DDT_Uspec[din][i] > 0])
            dout = choices(range(len(DDT_Uspec[din])), weights=DDT_Uspec[din], k=1)[0]
            douts_U.append(dout)
    
    dout_S = tip.compose(douts_U)
    if debug:
        print(f"{din_S} = {dins_U} -> {douts_U} = {dout_S}")
    return din_S, dout_S

def valid_din_dout_S(tip, din_S, dout_S, debug=False):
    """
    Checks if a given input-output difference pair `(din_S, dout_S)` is valid under the specified constraints.
    Note that multiplication with tip.R (Montgomery constant) is considered outside of this function.
    
    Additional constraints: 
    - For a tuple (xi,xi+di) with xi < xi+di, it shall hold that U(xi) < U(xi + di).

    Args:
    - tip (object): A member of the Tip5 family, used to decompose input and output differences.
    - din_S (int): The input difference for which validity is checked, din_S ~ [d1,d2,...,d8].
    - dout_S (int): The output difference associated with `din_S` to check for validity.
    - debug (bool, optional): If set to True, the function may print additional debug information. Defaults to False.

    Returns:
    - bool: `True` if the transition from `din_S` to `dout_S` is valid, `False` otherwise.
    """
    dins_U = tip.decompose(din_S) # [d1,d2,...,d8]
    douts_U = tip.decompose(dout_S)
    
    for din, dout in zip(dins_U, douts_U):
        if DDT_Uspec[din][dout] == 0:
            return False
        
    return True

def din_dout_to_in_S(tip, din_S, dout_S, debug=False):
    """
    Given a valid input-output difference pair `(din_S, dout_S)`, finds possible input values such that 
    the difference in the S-box output matches the given output difference. 
    Note that multiplication with tip.R (Montgomery constant) is applied to inputs to conform to original description of S-box.
    
    Additional constraints: 
    - For a tuple (xi,xi+di) with xi < xi+di, it shall hold that U(xi) < U(xi + di).
    
    Args:
    - tip (object): A member of the Tip5 family, used to decompose and compose differences.
    - din_S (int): The input difference to the S-box, din_S ~ [d1,d2,...,d8].
    - dout_S (int): The output difference to the S-box.
    - debug (bool, optional): If True, prints intermediate steps for debugging purposes. Defaults to False.
    
    Returns:
    - list: A list of possible input values (`ins_S`) corresponding to the input-output difference pair `(din_S, dout_S)`, 
            that is, `S(x+din_S) - S(x) = dout_S` for all x in `ins_S`.
    """
    dins_U = tip.decompose(din_S) # [d1,d2,...,d8]
    douts_U = tip.decompose(dout_S)
    
    if debug:
        print(f"Differences U: {dins_U} -> {douts_U}")
    
    ins_U = []
    for din, dout in zip(dins_U, douts_U):
        assert(DDT_Uspec[din][dout] > 0) # should have been chosen like that
        ins_U.append([])
        if din == 0:
            assert(dout == 0)
            ins_U[-1].append(randint(0,2**8-1))
        else:
            # Find possible xi such that U(xi+din) - U(xi) = dout
            for xi in range(2**8-din):
                if tip.lut(xi + din) - tip.lut(xi) == dout: # this also ensures U(xi) < U(xi + di)
                    assert(dout == TIP5_LOOKUP_TABLE[xi + din] - TIP5_LOOKUP_TABLE[xi])
                    ins_U[-1].append(xi)
            assert(len(ins_U[-1]) == DDT_Uspec[din][dout])
    
    # Combine (cartesian product of) parts to value(s) in F
    ins_S = [tip.compose(list(in_S)) for in_S in list(product(*ins_U))]

    if debug:
        for in_S in ins_S:
            print('-'*80)
            # Generate concrete input to S
            x1 = tip.R_inv * in_S
            x2 = tip.R_inv * (in_S + din_S)

            print("Concrete input to S-Box")
            print(f"in_S = {x1}")
            print(f"in_S + din_S = {x2}")
            print(f"Difference: {x2-x1} = {tip.R_inv * din_S}")

            # This is what S does
            x1 = tip.R * x1 
            x2 = tip.R * x2

            parts_x1 = tip.decompose(x1)
            parts_x2 = tip.decompose(x2)

            print(f"R * in_S = {x1} = {parts_x1}")
            print(f"R * (in_S + din_S) = {x2} = {parts_x2}")
            print(f"Input difference: {x2-x1} = {[a-b for a,b in zip(parts_x2,parts_x1)]}")

            parts_y1 = [tip.lut(xi) for xi in parts_x1]
            parts_y2 = [tip.lut(xi) for xi in parts_x2]

            y1 = tip.compose(parts_y1)
            y2 = tip.compose(parts_y2)

            print(f"L(R * in_S) = {y1} = {parts_y1}")
            print(f"L(R * (in_S + din_S)) = {y2} = {parts_y2}")
            print(f"Output difference: {y2-y1} = {[a-b for a,b in zip(parts_y2,parts_y1)]}")

            y1 = tip.R_inv * y1
            y2 = tip.R_inv * y2

            print("Concrete output to S-Box")
            print(f"S(in_S) = {y1}")
            print(f"S(in_S + din_S) = {y2}")
            print(f"Difference: {y2-y1} = {tip.R_inv * dout_S}")
    
    assert(all([tip.S(tip.R_inv*(in_S + din_S)) - tip.S(tip.R_inv*in_S) == tip.R_inv*dout_S for in_S in ins_S]))
    
    return ins_S # possible x values

def din_to_dout_T(tip, din_T, debug=False):
    """
    For a given input difference `din_T` to a T-Box, computes the valid output differences.

    Args:
    - tip (object): A member of the Tip5 family.
    - din_T (int): The input difference for which the corresponding output difference is computed.
    - debug (bool, optional): If set to True, the function may print additional debug information. Defaults to False.

    Returns:
    - (din_T, dout_T) tuple:
        - din_T (int): The input difference.
        - dout_T (int): A valid output difference.
    """
    
    # Only non-zero difference ok here
    if din_T == 0:
        return 0, 0
    
    x = tip.random_element()
    dout_T = tip.T(x + din_T) - tip.T(x)
    
    return din_T, dout_T

def din_dout_to_in_T(tip, din_T, dout_T, debug=False):
    """
    For a given input-output difference pair (`din_T`, `dout_T`), returns all possible input values `x` such that 
    the T-Box transformation satisfies `T(x + din_T) - T(x) = dout_T`.

    Args:
    - tip (object): A member of the Tip5 family.
    - din_T (int): The input difference for which valid input values are computed.
    - dout_T (int): The output difference for which the corresponding input values are sought.
    - debug (bool, optional): If set to True, the function may print additional debug information. Defaults to False.

    Returns:
    - list: A list of possible input values `x` that satisfy the equation `T(x + din_T) - T(x) = dout_T`. 
            If no valid solutions exist, returns an empty list.
    """
    P_uni = PolynomialRing(tip.field, 'x')
    x = P_uni.gen()
    f = tip.T(x + din_T) - tip.T(x) - dout_T
    
    if f == 0: # This happens for zero-difference. Use random value.
        return [tip.random_element()]
    
    roots = f.roots(multiplicities=False)
    if len(roots) == 0: # No solution here
        return []
    
    return roots