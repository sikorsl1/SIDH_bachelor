import Arithmetic as arithm
import numpy as np
import Velu_algorithm as velAlg
import Serialization_tools as serialTool

#returns isogeny from given rational functions
def isogeny(alpha,beta,P):
    return np.array([alpha(P[0],P[1]),beta(P[0],P[1])])

#public parameters
p=431
a1 = arithm.ZPIfield(0,0,p)
a2 = arithm.ZPIfield(0,0,p)
a3 = arithm.ZPIfield(0,0,p)
a4 = arithm.ZPIfield(182,0,p)
a6 = arithm.ZPIfield(321,0,p)
coeffs = np.array([a1,a2,a3,a4,a6])
E = arithm.EllipticCurvePoint(coeffs,p,0)
Pa = np.array([arithm.ZPIfield(13,279,p),arithm.ZPIfield(49,119,p)])
Qa = np.array([arithm.ZPIfield(26,375,p),arithm.ZPIfield(151,166,p)])
Pb = np.array([arithm.ZPIfield(0,152,p),arithm.ZPIfield(154,420,p)])
Qb = np.array([arithm.ZPIfield(428,241,p),arithm.ZPIfield(248,146,p)])
print('--------------------')
print('Public parameters')
print(f'Prime number: {str(p)}')
print(f'Elliptic curve: {E.equation()}')
print(f'j-invariant of this elliptic curve: {str(E.j_invariant())}')
print(f'Alice\'s generators: Pa = {str(Pa)}, Qa = {str(Qa)}')
print(f'Bob\'s generators: Pb = {str(Pb)}, Qb = {str(Qb)}')
print('--------------------')

#secret parameters
n_a = 2
m_a = 11
n_b = 16
m_b = 1
A = E.add(E.m_mltpl(n_a,Pa),E.m_mltpl(m_a,Qa))
B = E.add(E.m_mltpl(n_b,Pb),E.m_mltpl(m_b,Qb))

#isogenies from base elliptic curve
A_cyclic = E.generate_cyclic_subroup(A)
B_cyclic = E.generate_cyclic_subroup(B)
V_A = velAlg.Velu_set(A_cyclic,E)
V_B = velAlg.Velu_set(B_cyclic,E)
new_coeffs_A, aux_values_A = velAlg.Velu_values(V_A,E)
new_coeffs_B, aux_values_B = velAlg.Velu_values(V_B,E)
alpha_x = velAlg.Velu_x(aux_values_A,E)
alpha_y = velAlg.Velu_y(aux_values_A,E)
beta_x = velAlg.Velu_x(aux_values_B,E)
beta_y = velAlg.Velu_y(aux_values_B,E)

#new elliptic curves and images of Alice's and Bob's secrects
E_A = arithm.EllipticCurvePoint(new_coeffs_A,p,0)
E_B = arithm.EllipticCurvePoint(new_coeffs_B,p,0)
Pa_beta = isogeny(beta_x,beta_y,Pa)
Qa_beta = isogeny(beta_x,beta_y,Qa)
Pb_alpha = isogeny(alpha_x,alpha_y,Pb)
Qb_alpha = isogeny(alpha_x,alpha_y,Qb)
A_beta = E_B.add(E_B.m_mltpl(n_a,Pa_beta),E_B.m_mltpl(m_a,Qa_beta))
B_alpha = E_A.add(E_A.m_mltpl(n_b,Pb_alpha),E_A.m_mltpl(m_b,Qb_alpha))
print('Image curves')
print(f'Alice\'s image elliptic curve: {E_A.equation()}')
print(f'j-invariant of this elliptic curve: {str(E_A.j_invariant())}')
print(f'Bob\'s image elliptic curve: {E_B.equation()}')
print(f'j-invariant of this elliptic curve: {str(E_B.j_invariant())}')
print('--------------------')

A_beta_cyclic = E_B.generate_cyclic_subroup(A_beta)
B_alpha_cyclic = E_A.generate_cyclic_subroup(B_alpha)
V_A_beta = velAlg.Velu_set(A_beta_cyclic,E_B)
V_B_alpha = velAlg.Velu_set(B_alpha_cyclic,E_A)
final_coeffs_AB,aux_values_A_beta = velAlg.Velu_values(V_A_beta,E_B)
final_coeffs_BA,aux_values_B_alpha = velAlg.Velu_values(V_B_alpha,E_A)
E_AB = arithm.EllipticCurvePoint(final_coeffs_AB,p,0)
E_BA = arithm.EllipticCurvePoint(final_coeffs_BA,p,0)
print('Shared secret')
print(f'j-invariants of final Alice\'s elliptic curve: {str(E_AB.j_invariant())}')
print(f'j-invariants of final Bob\'s elliptic curve: {str(E_BA.j_invariant())}')
print('--------------------')

#Little example of how to serialize and deserialize sets of points
# points_from_file = serialTool.deserialize_points('points431.npy',p)
# print(len(points_from_file))
# serialTool.serialize_points('points431_test.npy',points_from_file)