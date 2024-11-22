'''
Author: Fukang Liu.

This is the script to compute the coefficient matrix after Gaussian elimination for Tip5 and Monolith31.

The obtained matrices will be used in the file: Diff_Property_of_Sbox_and_Finding_Diffs.cpp.

In more details, they will be used in the following two functions:

1. findDiffTip5(int threadNum) 

2. findDiffMonolith31(int threadNum)
'''

from scipy.linalg import circulant

def computeCoeffMatrix_Monolith31():
    p = 2**31 - 1
    F = GF(p)
    _MDS32 = [0x536C316, 0x1DD20A84, 0x43E26541, 0x52B22B8D, 0x37DABDF0, 0x540EC006, 0x3015718D, 0x5A99E14C,

          0x23637285, 0x4C8A2F76, 0x5DEC4E6E, 0x374EE8D6, 0x27EDA4D8, 0x665D30D3, 0x32E44597, 0x43C7E2B3,

          0x67C4C603, 0x78A8631F, 0x452F77E3, 0x39F03DF, 0x743DBFE0, 0x4DA05A48, 0x5F027940, 0x8293632,

          0x50F2C76A, 0x7B773729, 0x577DE8B0, 0x73B1EAC6, 0x58DA7D29, 0x67AA4375,0xDBA9E33, 0x2655E5A1]


    M = Matrix(GF(2**31 - 1), circulant(_MDS32)).transpose()

    MT = M.submatrix(row=0,col=0,nrows=24,ncols=24)
    MTInv = MT.inverse()
    
    xM = matrix(GF(p),8,9, lambda i,j: MTInv[i+16][j] if j<8 else MTInv[i+16][23])
    xM_ChangeCols  = matrix(GF(p),8,9, lambda i,j: xM[i][(j+1)%9])
    xM_ef=xM_ChangeCols.echelon_form()

    yM = matrix(GF(p),8,10, lambda i,j: MT[i][j] if j<9 else MT[i][23])
    yM_ChangeCols  = matrix(GF(p),8,10, lambda i,j: yM[i][(j+1)%10]) 
    yM_ef=yM_ChangeCols.echelon_form()
    
    xCoeff = matrix(GF(p),8,1, lambda i,j: xM_ef[i][j+8])
    yCoeff= matrix(GF(p),8,2, lambda i,j: yM_ef[i][j+8])
    
    print("Coefficient matrix (last column) after Gaussian elimination for eqs in x")
    print(xCoeff)
    print()
    
    print("Coefficient matrix (last two columns) after Gaussian elimination for eqs in y")
    print(yCoeff)
    print()
    
    '''
    #verification
    e0 = (0x4a9f18c5-0x69737ce4)%p
    sol5 = matrix(GF(p),24,1, lambda i,j: e0 if i<1 else (-(xCoeff[i-1][0]*e0) if i<8 else 0))
    sol5Full = matrix(GF(p),24,1, lambda i,j: sol5[i][0] if i<23 else -(xCoeff[7][0]*e0))
    print("Delta x:")
    print(sol5Full.transpose())
    print()
    print("M^{-1}(Delta x):")
    print(MTInv*sol5Full)
    print()
    
    #verification
    #e1 = F.random_element()
    #e2 = F.random_element()
    e1 = 0x25cab5b6
    e2 = 0x267acfcc
    #print(hex(e1),hex(e2))
    sol6 = matrix(GF(p),24,1, lambda i,j: e1 if i<1 else (0 if i<23 else e2))
    sol6Full = matrix(GF(p),24,1, lambda i,j: sol6[i][0] if i%23<1 else (-(yCoeff[i-1][0]*e2+yCoeff[i-1][1]*e1) if i<9 else 0))
    print("Delta y:")
    print(sol6Full.transpose())
    print()
    print("M(Delta y):")
    print(MT*sol6Full)
    '''
    
def computeCoeffMatrix_Tip5():
    p = 2**64 - 2**32 + 1
    F = GF(p)
    m = 16
    R = (2**64)%p
    RIn = R.inverse_mod(p)
    #print(R,RIn, (R*RIn).inverse_mod(p))
    
    RM = diagonal_matrix(GF(p),m,[R,R,R,R,1,1,1,1,1,1,1,1,1,1,1,1])
    RInM = diagonal_matrix(GF(p),m,[RIn,RIn,RIn,RIn,1,1,1,1,1,1,1,1,1,1,1,1])
    
    MDS_MATRIX_FIRST_COLUMN = [

    61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689,

    26798, 17845
    
    ]
    
    M = matrix(GF(p), nrows=m, ncols=m, entries=circulant(MDS_MATRIX_FIRST_COLUMN[:m]))
    
    M = RM*M*RInM
    MIn = M.inverse()
    
    xM = matrix(GF(p),4,5, lambda i,j: MIn[i+12][j])
    yM = matrix(GF(p),4,5, lambda i,j: M[i][j])
    
    xM_ef = xM.echelon_form()
    yM_ef = yM.echelon_form()

    #print(xM_coeff)
    #print(yM_coeff)
    
    xM_coeff = matrix(GF(p),4,1, lambda i,j: p - xM_ef[i][j+4])
    yM_coeff = matrix(GF(p),4,1, lambda i,j: p-yM_ef[i][j+4])
    
    print("Coefficient matrix (last column) after Gaussian elimination for eqs in x")
    print(xM_coeff)
    print()
    
    print("Coefficient matrix (last column) after Gaussian elimination for eqs in y")
    print(yM_coeff)
    print()
   
    '''
    #verification
    e0 = 0xca4c16b2172b8982
    e1 = 0xfb06f2832291a5d9
    x = matrix(GF(p),16,1, lambda i,j: xM_coeff[i][0]*e0 if i<4 else(e0 if i==4 else 0))
    y = matrix(GF(p),16,1, lambda i,j: yM_coeff[i][0]*e1 if i<4 else(e1 if i==4 else 0))
    print("output diff of inner part: M^{-1}(Delta x)")
    print(MIn*x)
    print()
    
    print("input diff of the 3rd round: M(Delta y)")
    print(M*y)
    '''
    
if __name__ == "__main__":
    print("Monolith31:")
    computeCoeffMatrix_Monolith31()
    
    print()
    print("Tip5")
    computeCoeffMatrix_Tip5()




