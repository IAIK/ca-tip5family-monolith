def gb_attack_sage(eqs):
    """
    Computes the Gröbner basis and variety of a given polynomial system.

    Args:
    - eqs (list): A list of polynomial equations forming the system to be solved. The equations 
                  should be elements of the same polynomial ring.

    Returns:
    - tuple:
        - gblex (list): The computed Gröbner basis of the input system in lexicographic order.
        - V (list or None): The variety of solutions for the system if computable, otherwise `None`.
          If a `ValueError` occurs during the variety computation, it will be printed and `V` will be `None`.
    """
    P_lex = eqs[0].parent().change_ring(order='lex')
    eqs_lex = [P_lex(f) for f in eqs]
    gblex = Ideal(eqs_lex).groebner_basis()
    
    V = None
    try:
        V = Ideal(gblex).variety()
    except ValueError as e:
        print("A ValueError occurred during Variety computation:", e)
        
    return gblex, V

def convert_polynomial(f, substitution_dict, retained_vars):
    """
    Converts a polynomial `f` by substituting variables according to `substitution_dict` 
    and retaining only the specified variables in `retained_vars`.

    Parameters:
    - f (Polynomial): The original polynomial to be transformed.
    - substitution_dict (dict): A dictionary mapping variables in `f` to their new values for substitution.
    - retained_vars (list of str): Names of variables to retain in the resulting polynomial.

    Returns:
    Polynomial: The transformed polynomial in the ring with only the retained variables.
    """
    # Perform the substitution
    f2 = f.substitute(substitution_dict)
    
    # Define the new polynomial ring with only the retained variables
    R2 = PolynomialRing(f.base_ring(), names=retained_vars)
    
    # Extract the retained generators
    retained_gens = [f.parent().gens().index(var) for var in retained_vars]
    
    # Construct a polynomial in the new ring R2 using only retained variables
    f2_in_R2 = R2(0)
    for monomial, coeff in f2.dict().items():
        # Filter out the terms that involve only the retained variables
        new_monomial = tuple(monomial[i] for i in retained_gens)
        f2_in_R2 += R2({new_monomial: coeff})
    
    return f2_in_R2

def gaussian_elimination(eqs, debug=True):
    """
    Performs Gaussian elimination on a system of polynomial equations to separate linear 
    and non-linear equations, then applies substitution for solutions from the linear system.

    Parameters:
    - eqs (list of Polynomial): A list of polynomial equations to process.

    Returns:
    - tuple:
        - E3 (Sequence): A sequence of transformed non-linear equations after substitution.
        - sub_lin (dict): A dictionary of substitutions derived from solving the linear equations.
    """
    F = eqs[0].base_ring()
    R = eqs[0].parent()
    
    E = Sequence(eqs)
    E1 = Sequence([f for f in E if f.degree() == 1])
    E2 = Sequence([f for f in E if f.degree() != 1])
    
    if debug:
        print(f"Total: {len(E)} equations in {E.nvariables()} variables", flush=True)
        print(f"{len(E1)} linear equations in {E1.nvariables()} variables")
        print(f"{len(E2)} non-linear equations in {E2.nvariables()} variables")
    
    # Gaussian elimination and substitution of solutions from linear system
    M,v = E1.coefficient_matrix()
    sub_lin = {}
    for i, f in enumerate(M.echelon_form() * v):
        sub_lin[v[i][0]] = (v[i] - f)[0]
    
    new_vars = R.gens()[i+1:]
    E3 = Sequence([convert_polynomial(f,sub_lin,new_vars) for f in E2])
    
    if debug:
        print(f"-> {len(E3)} non-linear equations in {E3.nvariables()} variables")
    
    return E3, sub_lin

def solve_linear_equation_system(system):
    """
    Solves a system of linear equations where each equation is of the form `var + constant`. 
    The function assumes the system is already reduced and that the polynomials are univariate and of degree 1.

    Args:
    - system (Sequence[Polynomial]): A sequence of polynomials representing the linear equations.
      Each polynomial should be in the form `var + constant` (i.e., univariate polynomials of degree 1).

    Returns:
    - tuple: A tuple containing:
        - A dictionary mapping variable names (as monomials) to their corresponding solutions.
        - The reduced system (a sequence of polynomials after the reduction process).
    
    Raises:
    - AssertionError: If the system is not square or if any polynomial in the system is not univariate or not of degree 1.
    """
    assert(len(system) == system.nvariables()) # System given as sequence
    solutions = {}
    system = system.reduced() # This yields polynomials of the form: var + value
    assert(all([f.is_univariate() and f.degree() == 1 for f in system]))
    for f in system:
        if f == 0: # Exclude this case
            return {}, system
        v = f.monomials()[0] 
        c = f.constant_coefficient()
        solutions[v] = -c
    #print(solutions)
    return solutions, system