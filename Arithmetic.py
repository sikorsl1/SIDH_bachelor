import numpy as np
import Serialization_tools as st

#implementation of arithmetic in Z_p field
class ZPfield():

    def __init__(self,x,p):
        self.p = p
        self.x = x%self.p

    def __add__(self, other):
        return ZPfield((self.x+other.x)%self.p,self.p)

    def __sub__(self, other):
        return ZPfield((self.x-other.x)%self.p,self.p)

    def __mul__(self, other):
        return ZPfield((self.x * other.x)%self.p,self.p)

    def __eq__(self, other):
        return True if self.x == other.x else False

    def __str__(self):
        return str(self.x)

    def __repr__(self):
        return str(self.x)

    def __hash__(self):
        return hash(self.x)

    def egcd(self,a, b):
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = self.egcd(b % a, a)
            return (g, x - (b // a) * y, y) # // oznacza dzielenie całkowitoliczbowe

    def inv(self):
        g, x, y = self.egcd(self.x, self.p)
        if g != 1:
            raise Exception('modular inverse does not exist')
        else:
            return ZPfield(x % self.p,self.p)

#implementation of arithmetic in Z_p[i] field
class ZPIfield():

    #element in Zp[i] field has the form 'a+bi'
    def __init__(self, a, b, p):
        self.p = p
        self.x = np.array([ZPfield(a,p),ZPfield(b,p)])

    def __add__(self, other):
        tmp = self.x + other.x
        return ZPIfield(tmp[0].x,tmp[1].x,self.p)

    def __sub__(self, other):
        tmp = self.x - other.x
        return ZPIfield(tmp[0].x,tmp[1].x,self.p)

    def __mul__(self, other):
        tmp_a = self.x[0]*other.x[0]-self.x[1]*other.x[1]
        tmp_b = self.x[0]*other.x[1]+self.x[1]*other.x[0]
        return ZPIfield(tmp_a.x,tmp_b.x,self.p)

    def __str__(self):
        return str(self.x)

    def __repr__(self):
        return str(self.x)

    def __eq__(self, other):
        return True if self.x[0].x==other.x[0].x and self.x[1].x==other.x[1].x else False

    def __hash__(self):
        return hash(str(self.x))

    def inv(self):
        square = self.x*self.x
        denom = (square[0]+square[1]).inv()
        tmp_a = self.x[0]*denom
        tmp_b = (ZPfield(0,self.p)-self.x[1])*denom
        return ZPIfield(tmp_a.x,tmp_b.x,self.p)

    def const(self,c):
        return ZPIfield(c,0,self.p)

    #returns set of all squares in the field
    def all_squares(self):
        elements = [ZPIfield(i, j, self.p) for i in range(self.p) for j in range(self.p)]
        elements_var = {}
        for e in elements:
            if e*e in elements_var:
                elements_var[e*e].append(e)
            else :
                elements_var[e*e] = [e]
        return elements_var

class EllipticCurvePoint():

    def __init__(self,coeffs,p,generate_squares):
        self.coeffs = coeffs
        self.p = p
        if generate_squares == 1:
            self.squares = self.coeffs[0].all_squares()

    def lambdaVar(self,P,Q):
        if P[0]==Q[0]:
            return (self.const(3)*P[0]*P[0]+self.const(2)*self.coeffs[1]*P[0]+self.coeffs[3]-self.coeffs[0]*P[1])*\
                   ((self.const(2)*P[1]+self.coeffs[0]*P[0]+self.coeffs[2]).inv())
        else:
            return (Q[1]-P[1])*((Q[0]-P[0]).inv())

    def muVar(self,P,Q):
        if P[0]==Q[0]:
            return (self.const(0)-P[0]*P[0]*P[0]+self.coeffs[3]*P[0]+self.const(2)*self.coeffs[4]-self.coeffs[2]*P[1])*\
                   ((self.const(2)*P[1]+self.coeffs[0]*P[0]+self.coeffs[2]).inv())
        else:
            return (P[1]*Q[0]-P[0]*Q[1])*((Q[0]-P[0]).inv())

    def const(self,c):
        return self.coeffs[0].const(c)

    def add(self,P,Q):
        lambda_var = self.lambdaVar(P,Q)
        x = lambda_var*lambda_var
        x = x + self.coeffs[0]*lambda_var-self.coeffs[1]-P[0]-Q[0]
        y = self.const(0)-(lambda_var+self.coeffs[0])*x - self.muVar(P,Q) - self.coeffs[2]
        return np.array([x,y])

    def opposite(self,P):
        x = P[0]
        y = self.const(0)-P[1]-self.coeffs[0]*P[0]-self.coeffs[2]
        return np.array([x,y])

    #evaluates value of Weierstrass equation for given point
    def eval_weierstrass(self,P):
        tmp = P[1]*P[1]+self.coeffs[0]*P[0]*P[1] + self.coeffs[2]*P[1]-P[0]*P[0]*P[0]-\
              self.coeffs[1]*P[0]*P[0]-self.coeffs[3]*P[0]-self.coeffs[4]
        return tmp

    #function returns all points on elliptic curve
    #it works only for elliptic curve given by normal form of Weierstrass equation
    def find_points(self):
        elements = [ZPIfield(i, j, self.p) for i in range(self.p) for j in range(self.p)]
        points = []
        num_elements = len(elements)
        for index,e in enumerate(elements):
            var = e*e*e + e*self.coeffs[3] + self.coeffs[4]
            if var in self.squares :
                for e_i in self.squares[var]:
                    points.append(np.array([e,e_i]))
            if index%100==0:
                print(str(index*100//num_elements) + '% elementow przetworzono')
        return np.array(points)

    def equation(self):
        return 'y^2 + (' + str(self.coeffs[0]) + ')xy + (' + str(self.coeffs[2]) +\
               ')y = x^3 + (' + str(self.coeffs[1]) +  ')x^2 + (' + str(self.coeffs[3]) + ')x + (' + str(self.coeffs[4]) + ')'

    #returns m-th multiplicity of given point
    def m_mltpl(self,m, P):
        Q = P
        minus_P = self.opposite(P)
        itr = iter(range(m - 1))
        for i in itr:
            if Q[0] == minus_P[0] and Q[1] == minus_P[1]:
                if i == m - 2:
                    return 0 #this means that function returns point at infinity
                Q = P
                next(itr, None)
            else:
                Q = self.add(P, Q)
        return Q

    #generate cyclic subgroup generated by given point
    def generate_cyclic_subroup(self, P):
        Q = P
        subgroup = []
        minus_P = self.opposite(P)
        while Q[0] != minus_P[0] or Q[1] != minus_P[1]:
            subgroup.append(Q)
            Q = self.add(Q, P)
        subgroup.append(Q)
        return np.array(subgroup)

    def secret_cyclic_subgroup(P, rank, path, el):
        cyclic_subgroup = [el.m_mltpl(i + 1, P) for i in range(rank-1)]
        st.serialize_points(path, np.array(cyclic_subgroup))

    #evaluates rank of given point
    def rank_of_the_point(self,R):
        Q = R
        minus_R = self.opposite(R)
        i = 2
        while Q[0] != minus_R[0] or Q[1] != minus_R[1]:
            i = i + 1
            Q = self.add(Q, R)
        print(i)
        return i

    #generates all points of given rank
    def points4rank(self, m, points):
        subgroup = []
        points_num = len(points)
        for index,P in enumerate(points):
            minus_P = self.opposite(P)
            Q = self.m_mltpl(m - 1, P)
            if Q != 0:
                if minus_P[0] == Q[0] and minus_P[1] == Q[1]:
                    subgroup.append(P)
            if index % 100 == 0:
                print('Przetworzono ' + str(100 * index / points_num) + '% punktow.')
        return np.array(subgroup)

    def j_invariant(self):
        b1 = self.coeffs[0]*self.coeffs[0] + self.const(4)*self.coeffs[1]
        b4 = self.coeffs[0]*self.coeffs[2] + self.const(2)*self.coeffs[3]
        b6 = self.coeffs[2]*self.coeffs[2] + self.const(4)*self.coeffs[4]
        l1 = (b1*b1 - self.const(24)*b4)*(b1*b1 - self.const(24)*b4)*(b1*b1 - self.const(24)*b4)
        nom = self.const(1728)*l1
        denom = l1 - (self.const(36)*b1*b4-b1*b1*b1-self.const(216)*b6)*(self.const(36)*b1*b4-b1*b1*b1-self.const(216)*b6)
        return nom*denom.inv()