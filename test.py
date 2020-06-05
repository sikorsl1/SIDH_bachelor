import Arithmetic as ar
import numpy as np
import Velu_algorithm as va

#returns isogeny from given rational functions
def isogeny(alpha,beta,P):
    return np.array([alpha(P[0],P[1]),beta(P[0],P[1])])

#public parameters
p=431
a1 = ar.ZPIfield(0,0,p)
a2 = ar.ZPIfield(0,0,p)
a3 = ar.ZPIfield(0,0,p)
a4 = ar.ZPIfield(182,0,p)
a6 = ar.ZPIfield(321,0,p)
coeffs = np.array([a1,a2,a3,a4,a6])
E = ar.EllipticCurvePoint(coeffs,p,0)
Pa = np.array([ar.ZPIfield(13,279,p),ar.ZPIfield(49,119,p)])
Qa = np.array([ar.ZPIfield(26,375,p),ar.ZPIfield(151,166,p)])
Pb = np.array([ar.ZPIfield(0,152,p),ar.ZPIfield(154,420,p)])
Qb = np.array([ar.ZPIfield(428,241,p),ar.ZPIfield(248,146,p)])
print('--------------------')
print('Public parameters')
print('Prime number: '+str(p))
print('Elliptic curve: ' + E.equation())
print('j-invariant of this elliptic curve: ' + str(E.j_invariant()))
print('Alice\'s generators: Pa = ' + str(Pa) + ', Qa = ' + str(Qa))
print('Bob\'s generators: Pb = ' + str(Pb) + ', Qb = ' + str(Qb))
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
V_A = va.Velu_set(A_cyclic,E)
V_B = va.Velu_set(B_cyclic,E)
new_coeffs_A, aux_values_A = va.Velu_values(V_A,E)
new_coeffs_B, aux_values_B = va.Velu_values(V_B,E)
alpha_x = va.Velu_x(aux_values_A,E)
alpha_y = va.Velu_y(aux_values_A,E)
beta_x = va.Velu_x(aux_values_B,E)
beta_y = va.Velu_y(aux_values_B,E)

#new elliptic curves and images of Alice's and Bob's secrects
E_A = ar.EllipticCurvePoint(new_coeffs_A,p,0)
E_B = ar.EllipticCurvePoint(new_coeffs_B,p,0)
Pa_beta = isogeny(beta_x,beta_y,Pa)
Qa_beta = isogeny(beta_x,beta_y,Qa)
Pb_alpha = isogeny(alpha_x,alpha_y,Pb)
Qb_alpha = isogeny(alpha_x,alpha_y,Qb)
A_beta = E_B.add(E_B.m_mltpl(n_a,Pa_beta),E_B.m_mltpl(m_a,Qa_beta))
B_alpha = E_A.add(E_A.m_mltpl(n_b,Pb_alpha),E_A.m_mltpl(m_b,Qb_alpha))
print('Image curves')
print('Alice\'s image elliptic curve: ' + E_A.equation())
print('j-invariant of this elliptic curve: ' + str(E_A.j_invariant()))
print('Bob\'s image elliptic curve: ' + E_B.equation())
print('j-invariant of this elliptic curve: ' + str(E_B.j_invariant()))
print('--------------------')

A_beta_cyclic = E_B.generate_cyclic_subroup(A_beta)
B_alpha_cyclic = E_A.generate_cyclic_subroup(B_alpha)
V_A_beta = va.Velu_set(A_beta_cyclic,E_B)
V_B_alpha = va.Velu_set(B_alpha_cyclic,E_A)
final_coeffs_AB,aux_values_A_beta = va.Velu_values(V_A_beta,E_B)
final_coeffs_BA,aux_values_B_alpha = va.Velu_values(V_B_alpha,E_A)
E_AB = ar.EllipticCurvePoint(final_coeffs_AB,p,0)
E_BA = ar.EllipticCurvePoint(final_coeffs_BA,p,0)
print('Shared secret')
print('j-invariants of final Alice\'s elliptic curve: ' + str(E_AB.j_invariant()))
print('j-invariants of final Bob\'s elliptic curve: ' + str(E_BA.j_invariant()))
print('--------------------')