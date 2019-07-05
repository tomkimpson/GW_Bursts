MODULE init_EQL

  
  
USE SHARED_CONSTANTS
IMPLICIT NONE


PUBLIC :: EQL_circular, EQL_eccentric
PRIVATE :: func_f, func_g, func_h, func_d

CONTAINS

SUBROUTINE EQL_circular()
IMPLICIT NONE
real(kind=dp) :: N,dL



print *, '~ circular orbit with r = ', r_init




N = (1.00_dp - 3.00_dp/r_init + 2.00_dp*a*r_init**(-1.50_dp))**0.50_dp
E = (1.00_dp - 2.00_dp/r_init + a*r_init**(-1.50_dp))/N
L = r_init**0.50_dp * (1+(a/r_init)**2.00_dp - 2*a*r_init**(-1.50_dp))/N
dL = 0.00_dp
L = L +dL


Q = 0.00_dp



END SUBROUTINE EQL_circular


SUBROUTINE EQL_eccentric(rp,ra,theta)
real(kind=dp) rp,ra,theta
real(kind=dp) z, f1,g1,h1,d1,f2,g2,h2,d2
real(kind=dp) kappa,eps,rho,et,sig
real(kind=dp) E_top, E_bot,DD


if (a .LT. 0) then
DD = -1.0_dp
else
DD = +1.0_dp
endif








z = cos(theta)**2



print *, ra, rp, z


call func_f(rp,a,z,f1)
call func_g(rp,a,z,g1)
call func_h(rp,a,z,h1)
call func_d(rp,a,z,d1)

call func_f(ra,a,z,f2)
call func_g(ra,a,z,g2)
call func_h(ra,a,z,h2)
call func_d(ra,a,z,d2)

kappa = d1*h2 - d2*h1
eps = d1*g2 - d2*g1
rho = f1*h2 - f2*h1
et = f1*g2 - f2*g1
sig = g1*h2 - g2*h1
                        
E_top = kappa*rho + 2.0_dp*eps*sig - 2.0_dp*(sig*(sig*eps**2 + rho*eps*kappa-et*kappa**2))**0.5
E_bot = rho**2 + 4.0_dp*et*sig

E = (E_top/E_bot)**0.5
L = -g1*E/h1 + DD*(g1**2*E**2+(f1*E**2-d1)*h1)**0.5/h1
Q = z*(a**2*(1-E**2)+L**2/(1.0_dp-z))


print *, 'E,L,Q = ', E, L, Q
print *, 'Inclination angle i = ', atan(sqrt(Q)/L) * 180.0_dp/PI

stop


END SUBROUTINE EQL_eccentric



SUBROUTINE func_f(r,a,z,out)
real(kind=dp) r, a, z ,out
real(kind=dp) delta
delta = r**2.0 - 2.0*r +a**2.0
out = r**4 +a**2*(r*(r+2)+z*delta)
END SUBROUTINE func_f

SUBROUTINE func_g(r,a,z,out)
real(kind=dp) r, a, z, out
out = 2*a*r
END SUBROUTINE func_g


SUBROUTINE func_h(r,a,z,out)
real(kind=dp) r, a, z, out,delta
delta = r**2.0 - 2.0*r +a**2.0
out = r*(r-2)+z*delta/(1-z)
END SUBROUTINE func_h

SUBROUTINE func_d(r,a,z,out)
real(kind=dp) r,a,z,out,delta
delta = r**2.0 - 2.0*r +a**2.0
out = (r**2+a**2*z)*delta
END SUBROUTINE func_d












END MODULE init_EQL
