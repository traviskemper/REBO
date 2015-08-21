c-----------------------------------------------------------------------
c     bicub performs a bicubic spline interpolation.  two variables are
c     passed in -- x and y.  these are used with the array of 16
c     coefficients coeff to calculate the value of the spline f, as well
c     as the derivatives df/dx and df/dy.  if the values of x and y need
c     to be screened, this should be done BEFORE calling bicub -- this
c     routine knows nothing about the proper domain of the spline.
c     the only values changed in this routine are f, dfdx, and dfdy.
c     x and y are untouched.  this routine performs the same function as
c     the old routine bcuint() used to.
c     ...sjs 4/30/98
c-----------------------------------------------------------------------
c
      subroutine bicub(x, y, coeff, f, dfdx, dfdy)

      IMPLICIT none

      REAL*8 :: coeff(16)
     &,x,y,f,dfdx,dfdy
     &,temp,tempf,tempfy,tempfx

c     zero the output variables:

      f = 0.d0
      dfdx = 0.d0
      dfdy = 0.d0

c     f is generated from a 16-term polynomial.  dfdx and dfdy are
c     derivatives.  the form may look ugly, but it's fast

      temp = coeff(13) + x * (coeff(14) 
     .     + x * (coeff(15) + x * coeff(16)))
      tempf = temp
      tempfy = 3.d0 * temp

      temp = coeff(9) + x * (coeff(10)
     .     + x * (coeff(11) + x * coeff(12)))
      tempf = tempf * y + temp
      tempfy = tempfy * y + 2.d0 * temp

      temp = coeff(5) + x * (coeff(6)
     .     + x * (coeff(7) + x * coeff(8)))
      tempf = tempf * y + temp
      tempfy = tempfy * y + temp

      temp = coeff(1) + x * (coeff(2)
     .     + x * (coeff(3) + x * coeff(4)))
      tempf = tempf * y + temp

      f = tempf

      dfdy = tempfy

      temp = coeff(4) + y * (coeff(8)
     .     + y * (coeff(12) + y * coeff(16)))
      tempfx = 3.d0 * temp

      temp = coeff(3) + y * (coeff(7)
     .     + y * (coeff(11) + y * coeff(15)))
      tempfx = tempfx * x + 2.d0 * temp

      temp = coeff(2) + y * (coeff(6)
     .     + y * (coeff(10) + y * coeff(14)))
      tempfx = tempfx * x + temp

      dfdx = tempfx

      return
      end
c-----------------------------------------------------------------------
c     bcucof calculates the coefficients of a bicubic spline, given the
c     value of the spline and three of its derivatives at each of the
c     four knots.  Taken from Numerical Recipes....sjs 6/11/01
c-----------------------------------------------------------------------
c
      subroutine bcucof(xmin, xmax, ymin, ymax, y, y1, y2, y12, dl)

      implicit real*8 (a-h,o-z)

c     xmin, xmax, ymin, ymax are the knots
c     y(4) are the spline values at the knots
c     y1(4) are the d/dx values at the knot
c     y2(4) are the d/dy values at the knots
c     y12(4) are the d2/dxdy values at the knots
c     dl(16) are the spline coefficients

      INTEGER :: xmin, xmax, ymin, ymax

      REAL*8 :: c(4,4), c2(4,4), y(4), y1(4), y2(4), y12(4), 
     &     cl(16), dl(16), x(16)
      INTEGER ::  wt(16,16)

      data wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,
     .     10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,
     .     4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,
     .     10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,
     .     0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,
     .     10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,
     .     5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,
     .     10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

      dx = xmax - xmin
      dy = ymax - ymin

      dxdy = dx * dy

      do 11 i = 1, 4
         x(i) = y(i)
         x(i+4) = y1(i) * dx
         x(i+8) = y2(i) * dy
         x(i+12) = y12(i) * dxdy
 11   continue

      do 13 i = 1, 16
         xx = 0.d0
         do 12 k = 1, 16
            xx = xx + wt(i,k) * x(k)
 12      continue
         cl(i) = xx
 13   continue
      
      l = 0
      do 15 i = 1, 4
         do 14 j = 1, 4
            l = l + 1
            c(i,j) = cl(l)
 14      continue
 15   continue

      c2(1,1) = c(1,1)
c
      c2(1,2) = c(1,2) / dy
      c2(1,1) = c2(1,1) - c(1,2) * ymin / dy
c
      temp = c(1,3) / dy**2
      c2(1,3) = temp
      c2(1,2) = c2(1,2) - 2 * temp * ymin
      c2(1,1) = c2(1,1) + temp * ymin**2
c
      temp = c(1,4) / dy**3
      c2(1,4) = temp
      c2(1,3) = c2(1,3) - 3 * temp * ymin
      c2(1,2) = c2(1,2) + 3 * temp * ymin**2
      c2(1,1) = c2(1,1) - temp * ymin**3
c
      temp = c(2,1) / dx
      c2(2,1) = temp
      c2(1,1) = c2(1,1) - temp * xmin
c
      temp = c(2,2) / dx / dy
      c2(2,2) = temp
      c2(2,1) = c2(2,1) - temp * ymin
      c2(1,2) = c2(1,2) - temp * xmin
      c2(1,1) = c2(1,1) + temp * xmin * ymin
c
      temp = c(2,3) / dx / dy**2
      c2(2,3) = temp
      c2(2,2) = c2(2,2) - 2 * temp * ymin
      c2(2,1) = c2(2,1) + temp * ymin**2
      c2(1,3) = c2(1,3) - temp * xmin
      c2(1,2) = c2(1,2) + 2 * temp * xmin * ymin
      c2(1,1) = c2(1,1) - temp * xmin * ymin**2
c
      temp = c(2,4) / dx / dy**3
      c2(2,4) = temp
      c2(2,3) = c2(2,3) - 3 * temp * ymin
      c2(2,2) = c2(2,2) + 3 * temp * ymin**2
      c2(2,1) = c2(2,1) - temp * ymin**3
      c2(1,4) = c2(1,4) - temp * xmin
      c2(1,3) = c2(1,3) + 3 * temp * xmin * ymin
      c2(1,2) = c2(1,2) - 3 * temp * xmin * ymin**2
      c2(1,1) = c2(1,1) + temp * xmin * ymin**3
c
      temp = c(3,1) / dx**2
      c2(3,1) = temp
      c2(2,1) = c2(2,1) - 2 * temp * xmin
      c2(1,1) = c2(1,1) + temp * xmin**2
c
      temp = c(3,2) / dx**2 / dy
      c2(3,2) = temp
      c2(3,1) = c2(3,1) - temp * ymin
      c2(2,2) = c2(2,2) - 2 * temp * xmin
      c2(2,1) = c2(2,1) + 2 * temp * xmin * ymin
      c2(1,2) = c2(1,2) + temp * xmin**2
      c2(1,1) = c2(1,1) - temp * xmin**2 * ymin
c
      temp = c(3,3) / dx**2 / dy**2
      c2(3,3) = temp
      c2(3,2) = c2(3,2) - 2 * temp * ymin
      c2(3,1) = c2(3,1) + temp * ymin**2
      c2(2,3) = c2(2,3) - 2 * temp * xmin
      c2(2,2) = c2(2,2) + 4 * temp * xmin * ymin
      c2(2,1) = c2(2,1) - 2 * temp * xmin * ymin**2
      c2(1,3) = c2(1,3) + temp * xmin**2
      c2(1,2) = c2(1,2) - 2 * temp * xmin**2 * ymin
      c2(1,1) = c2(1,1) + temp * xmin**2 * ymin**2
c
      temp = c(3,4) / dx**2 / dy**3
      c2(3,4) = temp
      c2(3,3) = c2(3,3) - 3 * temp * ymin
      c2(3,2) = c2(3,2) + 3 * temp * ymin**2
      c2(3,1) = c2(3,1) - temp * ymin**3
      c2(2,4) = c2(2,4) - 2 * temp * xmin
      c2(2,3) = c2(2,3) + 6 * temp * xmin * ymin
      c2(2,2) = c2(2,2) - 6 * temp * xmin * ymin**2
      c2(2,1) = c2(2,1) + 2 * temp * xmin * ymin**3
      c2(1,4) = c2(1,4) + temp * xmin**2
      c2(1,3) = c2(1,3) - 3 * temp * xmin**2 * ymin
      c2(1,2) = c2(1,2) + 3 * temp * xmin**2 * ymin**2
      c2(1,1) = c2(1,1) - temp * xmin**2 * ymin**3
c
      temp = c(4,1) / dx**3
      c2(4,1) = temp
      c2(3,1) = c2(3,1) - 3 * temp * xmin
      c2(2,1) = c2(2,1) + 3 * temp * xmin**2
      c2(1,1) = c2(1,1) - temp * xmin**3
c
      temp = c(4,2) / dx**3 / dy
      c2(4,2) = temp
      c2(4,1) = c2(4,1) - temp * ymin
      c2(3,2) = c2(3,2) - 3 * temp * xmin
      c2(3,1) = c2(3,1) + 3 * temp * xmin * ymin
      c2(2,2) = c2(2,2) + 3 * temp * xmin**2
      c2(2,1) = c2(2,1) - 3 * temp * xmin**2 * ymin
      c2(1,2) = c2(1,2) - temp * xmin**3
      c2(1,1) = c2(1,1) + temp * xmin**3 * ymin
c
      temp = c(4,3) / dx**3 / dy**2
      c2(4,3) = temp
      c2(4,2) = c2(4,2) - 2 * temp * ymin
      c2(4,1) = c2(4,1) + temp * ymin**2
      c2(3,3) = c2(3,3) - 3 * temp * xmin
      c2(3,2) = c2(3,2) + 6 * temp * xmin * ymin
      c2(3,1) = c2(3,1) - 3 * temp * xmin * ymin**2
      c2(2,3) = c2(2,3) + 3 * temp * xmin**2
      c2(2,2) = c2(2,2) - 6 * temp * xmin**2 * ymin
      c2(2,1) = c2(2,1) + 3 * temp * xmin**2 * ymin**2
      c2(1,3) = c2(1,3) - temp * xmin**3
      c2(1,2) = c2(1,2) + 2 * temp * xmin**3 * ymin
      c2(1,1) = c2(1,1) - temp * xmin**3 * ymin**2
c
      temp = c(4,4) / dx**3 / dy**3
      c2(4,4) = temp
      c2(4,3) = c2(4,3) - 3 * temp * ymin
      c2(4,2) = c2(4,2) + 3 * temp * ymin**2
      c2(4,1) = c2(4,1) - temp * ymin**3
      c2(3,4) = c2(3,4) - 3 * temp * xmin
      c2(3,3) = c2(3,3) + 9 * temp * xmin * ymin
      c2(3,2) = c2(3,2) - 9 * temp * xmin * ymin**2
      c2(3,1) = c2(3,1) + 3 * temp * xmin * ymin**3
      c2(2,4) = c2(2,4) + 3 * temp * xmin**2
      c2(2,3) = c2(2,3) - 9 * temp * xmin**2 * ymin
      c2(2,2) = c2(2,2) + 9 * temp * xmin**2 * ymin**2
      c2(2,1) = c2(2,1) - 3 * temp * xmin**2 * ymin**3
      c2(1,4) = c2(1,4) - temp * xmin**3
      c2(1,3) = c2(1,3) + 3 * temp * xmin**3 * ymin
      c2(1,2) = c2(1,2) - 3 * temp * xmin**3 * ymin**2
      c2(1,1) = c2(1,1) + temp * xmin**3 * ymin**3

      k=0
      do 17 i = 1,4
         do 16 j = 1,4
            k = k + 1
            dl(k) = c2(i,j)
 16      continue
 17   continue

      end
c
c-----------------------------------------------------------------------
c     tricub performs a tricubic spline interpolation.  three variables
c     are passed in -- x, y, and z.  these are used with the array of
c     64 coefficients coeff to calculate the
c     value of the spline f, as well as the derivatives df/dx, df/dy,
c     and df/dz.
c     if the values of x, y, and z need to be screened, this should be
c     done BEFORE calling tricub -- this routine 
c     knows nothing about the proper domain of the spline.
c     the only values changed in this routine are f, dfdx, dfdy, and
c     dfdz.  x, y, and z are untouched.
c     this routine performs the same function as the old routines
c     tor() and radic() used to, although some of the screening done
c     by those routines is now done in pibond and/or spoof.
c     ...sjs 3/24/98
c-----------------------------------------------------------------------
c     A completely general version of this routine should perhaps
c     return d^2f/dxdy, d^2f/dxdz, d^2f/dydz, and d^3f/dxdydz, since
c     these can be used in constructing the coefficients.  But nobody
c     needs these values currently so they are not returned.
c-----------------------------------------------------------------------
c     
      subroutine tricub(x, y, z, coeff, f, dfdx, dfdy, dfdz)

      implicit real*8 (a-h,o-z)

      dimension coeff(64)

c     zero the output variables:

      f = 0.d0
      dfdx = 0.d0
      dfdy = 0.d0
      dfdz = 0.d0

c     f is generated from a 64-term polynomial.  dfdx, dfdy, and dfdz
c     are derivatives.  the form may look ugly, but it's fast

      temp = coeff(61) + z * (coeff(62) 
     .     + z * (coeff(63) + z * coeff(64)))
      temp2 = temp
      temp5 = temp
      temp = coeff(57) + z * (coeff(58)
     .     + z * (coeff(59) + z * coeff(60)))
      temp2 = temp2 * y + temp
      temp6 = temp
      temp = coeff(53) + z * (coeff(54)
     .     + z * (coeff(55) + z * coeff(56)))
      temp2 = temp2 * y + temp
      temp7 = temp
      temp = coeff(49) + z * (coeff(50)
     .     + z * (coeff(51) + z * coeff(52)))
      temp2 = temp2 * y + temp
      temp3 = temp2
      temp4 = 3.d0 * temp2

      temp = coeff(45) + z * (coeff(46)
     .     + z * (coeff(47) + z * coeff(48)))
      temp2 = temp
      temp5 = temp5 * x + temp
      temp = coeff(41) + z * (coeff(42)
     .     + z * (coeff(43) + z * coeff(44)))
      temp2 = temp2 * y + temp
      temp6 = temp6 * x + temp
      temp = coeff(37) + z * (coeff(38)
     .     + z * (coeff(39) + z * coeff(40)))
      temp2 = temp2 * y + temp
      temp7 = temp7 * x + temp
      temp = coeff(33) + z * (coeff(34)
     .     + z * (coeff(35) + z * coeff(36)))
      temp2 = temp2 * y + temp
      temp3 = temp3 * x + temp2
      temp4 = temp4 * x + 2.d0 * temp2

      temp = coeff(29) + z * (coeff(30)
     .     + z * (coeff(31) + z * coeff(32)))
      temp2 = temp
      temp5 = temp5 * x + temp
      temp = coeff(25) + z * (coeff(26)
     .     + z * (coeff(27) + z * coeff(28)))
      temp2 = temp2 * y + temp
      temp6 = temp6 * x + temp
      temp = coeff(21) + z * (coeff(22)
     .     + z * (coeff(23) + z * coeff(24)))
      temp2 = temp2 * y + temp
      temp7 = temp7 * x + temp
      temp = coeff(17) + z * (coeff(18)
     .     + z * (coeff(19) + z * coeff(20)))
      temp2 = temp2 * y + temp
      temp3 = temp3 * x + temp2
      temp4 = temp4 * x + temp2
      dfdx = dfdx + temp4

      temp = coeff(13) + z * (coeff(14)
     .     + z * (coeff(15) + z * coeff(16)))
      temp2 = temp
      temp5 = temp5 * x + temp
      temp = coeff(9) + z * (coeff(10)
     .     + z * (coeff(11) + z * coeff(12)))
      temp2 = temp2 * y + temp
      temp6 = temp6 * x + temp
      temp = coeff(5) + z * (coeff(6)
     .     + z * (coeff(7) + z * coeff(8)))
      temp2 = temp2 * y + temp
      temp7 = temp7 * x + temp
      temp = coeff(1) + z * (coeff(2)
     .     + z * (coeff(3) + z * coeff(4)))
      temp2 = temp2 * y + temp
      temp3 = temp3 * x + temp2
      f = f + temp3
      temp8 = temp7 + y * (2.d0 * temp6 + y * 3.d0 * temp5)
      dfdy = dfdy + temp8

      temp = coeff(52) + y * (coeff(56) 
     .     + y * (coeff(60) + y * coeff(64)))
      temp2 = temp
      temp = coeff(36) + y * (coeff(40) 
     .     + y * (coeff(44) + y * coeff(48)))
      temp2 = temp2 * x + temp
      temp = coeff(20) + y * (coeff(24) 
     .     + y * (coeff(28) + y * coeff(32)))
      temp2 = temp2 * x + temp
      temp = coeff(4) + y * (coeff(8) 
     .     + y * (coeff(12) + y * coeff(16)))
      temp2 = temp2 * x + temp
      temp3 = 3.d0 * temp2

      temp = coeff(51) + y * (coeff(55) 
     .     + y * (coeff(59) + y * coeff(63)))
      temp2 = temp
      temp = coeff(35) + y * (coeff(39) 
     .     + y * (coeff(43) + y * coeff(47)))
      temp2 = temp2 * x + temp
      temp = coeff(19) + y * (coeff(23) 
     .     + y * (coeff(27) + y * coeff(31)))
      temp2 = temp2 * x + temp
      temp = coeff(3) + y * (coeff(7) 
     .     + y * (coeff(11) + y * coeff(15)))
      temp2 = temp2 * x + temp
      temp3 = temp3 * z + 2.d0 * temp2

      temp = coeff(50) + y * (coeff(54) 
     .     + y * (coeff(58) + y * coeff(62)))
      temp2 = temp
      temp = coeff(34) + y * (coeff(38) 
     .     + y * (coeff(42) + y * coeff(46)))
      temp2 = temp2 * x + temp
      temp = coeff(18) + y * (coeff(22) 
     .     + y * (coeff(26) + y * coeff(30)))
      temp2 = temp2 * x + temp
      temp = coeff(2) + y * (coeff(6) 
     .     + y * (coeff(10) + y * coeff(14)))
      temp2 = temp2 * x + temp
      temp3 = temp3 * z + temp2
      dfdz = dfdz + temp3

      return 
      end
c
c-----------------------------------------------------------------------
c     tcucof calculates the coefficients of a tricubic spline, given
c     the values of the function (f), the first deriviatives in each 
c     direction (fx, fy, fz), the mixed second derivatives (fxy, fxz,
c     fyz), and the mixed third derivative (fxyz) at each of the
c     8 knots.
c     64x64 matrix used to obtain the coefficients on the unit cube
c     obtained by Yang Li.
c     ...sjs 2/11/02
c-----------------------------------------------------------------------
c
      subroutine tcucofair(xmin, xmax, ymin, ymax, zmin, zmax, 
     .     f, fx, fy, fz, fxy, fxz, fyz, fxyz, coeff2)

      implicit real*8 (a-h,o-z)

c     xmin, xmax, ymin, ymax, zmin, zmax are the knot coordinates
c     f is the function value at each knot
c     fx, fy, fz are the first derivatives at each knot
c     fxy, fxz, fyz are the 2nd derivatives at each knot
c     fxyz is the third derivative at each knot
c     coeff holds the spline coefficients

      dimension finput(0:1,0:1,0:1,1:8),
     .     coeff(0:3,0:3,0:3), coeff2(0:3,0:3,0:3),
     .     f(0:1,0:1,0:1), fx(0:1,0:1,0:1), fy(0:1,0:1,0:1), 
     .     fz(0:1,0:1,0:1), fxy(0:1,0:1,0:1), fxz(0:1,0:1,0:1),
     .     fyz(0:1,0:1,0:1), fxyz(0:1,0:1,0:1),
     .     onefac(0:9),
     .     dxfac(0:3), dyfac(0:3), dzfac(0:3),
     .     xmfac(0:3), ymfac(0:3), zmfac(0:3)

      dimension imat(0:1,0:1,0:1,1:8,0:3,0:3,0:3),
     .     ichoos(0:3,0:3)

c     ichoos(i,j) = "i choose j" = binomial coefficient
c        = i! / j! / (i-j)!
c     should only be evaluated for i >= j.

      data ichoos/1, 1, 1, 1, -1, 1, 2, 3, -1, -1, 1, 3, -1, -1 , -1, 1/

c     Mathematical summary:  there is a tricubic function,
c     f(x,y,z) = \sum_{i=0}^3 \sum_{j=0}^3 \sum_{k=0}^3
c                    c_{ijk} x^i y^j z^k
c     whose coefficients a_{ijk} are uniquely determined by specifying
c     the values of f(x,y,z) at the 8 knot points obtained by taking
c     every combination of {x_{min},x_{max}} with {y_{min},y_{max}} and
c     {z_{min},z_{max}}.  This routine calculates the c_{ijk}.

c     conceptually, imat is a 64x64 matrix that, when multiplied by
c     a length-64 vector of input values (8 fxn vals and derivatives
c     at each of the 8 knots; 64 = 8x8) gives the length-64 vector of
c     coefficients in the cubic spline (coefficients for 0th through
c     3rd power of each of 3 input variables; 64=4^3).

c     computationally, imat is stored as a 2x2x2x8x4x4x4 tensor that,
c     when multiplied by a 2x2x2x8 tensor of input values gives a
c     4x4x4 tensor of coefficients:
c     imat(k,j,i,d,n,m,l) = coefficient of f^(d)(i,j,k) used to
c     calculate coeff2(n,m,l) = c_{lmn}.

c     finput stores the input function values and derivatives, in a
c     packed format:
c     finput(k,j,i,d) = f^(d)(i,j,k), where
c     f^(1) = f; f^(2) = df/dx; f^(3) = df/dy; f^(4) = df/dz;
c     f^(5) = d^2f/dxdy; f^(6) = d^2f/dxdz; f^(7) = d^2f/dydz;
c     f^(8) = d^3f/dxdydz;
c     i = {0,1} is the {lower,upper} value of the x coordinate;
c     j = {0,1} is the {lower,upper} value of the y coordinate;
c     k = {0,1} is the {lower,upper} value of the z coordinate.

c     the coeff* arrays store the tricubic spline coefficients,
c     at various intermediate and final stages.  they are stored with
c     indices in backwards order because of FORTRAN's array
c     storage order:
c     coeff(n,m,l) = a_{lmn} or b_{lmn}
c     coeff2(n,m,l) = c_{lmn}

c     the FORTRAN 77 standard says no more than 19 continuation
c     lines... luckily the compilers appear to be much more generous

      data imat/1,87*0,1,39*0,-3,3,22*0,-2,-1,38*0,2,-2,
     .22*0,2*1,54*0,1,95*0,1,31*0,-3,3,30*0,-2,-1,30*0,2,
     .-2,30*0,2*1,14*0,-3,0,3,13*0,-2,0,-1,69*0,-3,0,3,21*0,
     .-2,0,-1,13*0,9,2*-9,9,12*0,6,-6,3,-3,4*0,6,3,-6,-3,
     .20*0,4,2*2,1,12*0,-6,2*6,-6,12*0,-4,4,-2,2,4*0,2*-3,
     .2*3,20*0,2*-2,2*-1,12*0,2,0,-2,13*0,1,0,1,69*0,2,0,
     .-2,21*0,1,0,1,13*0,-6,2*6,-6,12*0,-3,3,-3,3,4*0,-4,
     .-2,4,2,20*0,-2,-1,-2,-1,12*0,4,2*-4,4,12*0,2,-2,2,
     .-2,4*0,2*2,2*-2,20*0,4*1,20*0,1,95*0,1,31*0,-3,3,30*0,
     .-2,-1,30*0,2,-2,30*0,2*1,54*0,1,87*0,1,39*0,-3,3,22*0,
     .-2,-1,38*0,2,-2,22*0,2*1,14*0,-3,0,3,21*0,-2,0,-1,
     .69*0,-3,0,3,13*0,-2,0,-1,13*0,9,2*-9,9,20*0,6,-6,3,
     .-3,4*0,6,3,-6,-3,12*0,4,2*2,1,12*0,-6,2*6,-6,20*0,
     .-4,4,-2,2,4*0,2*-3,2*3,12*0,2*-2,2*-1,12*0,2,0,-2,
     .21*0,1,0,1,69*0,2,0,-2,13*0,1,0,1,13*0,-6,2*6,-6,20*0,
     .-3,3,-3,3,4*0,-4,-2,4,2,12*0,-2,-1,-2,-1,12*0,4,2*-4,
     .4,20*0,2,-2,2,-2,4*0,2*2,2*-2,12*0,4*1,4*0,-3,3*0,
     .3,3*0,-2,3*0,-1,75*0,-3,3*0,3,11*0,-2,3*0,-1,19*0,
     .9,-9,2*0,-9,9,2*0,6,-6,2*0,3,-3,10*0,6,3,2*0,-6,-3,
     .10*0,4,2,2*0,2,1,18*0,-6,6,2*0,6,-6,2*0,-4,4,2*0,-2,
     .2,10*0,2*-3,2*0,2*3,10*0,2*-2,2*0,2*-1,34*0,-3,3*0,
     .3,11*0,-2,3*0,-1,75*0,-3,3*0,3,3*0,-2,3*0,-1,19*0,
     .9,-9,2*0,-9,9,10*0,6,-6,2*0,3,-3,10*0,6,3,2*0,-6,-3,
     .2*0,4,2,2*0,2,1,18*0,-6,6,2*0,6,-6,10*0,-4,4,2*0,-2,
     .2,10*0,2*-3,2*0,2*3,2*0,2*-2,2*0,2*-1,2*0,9,0,-9,0,
     .-9,0,9,0,6,0,-6,0,3,0,-3,0,6,0,3,0,-6,0,-3,9*0,4,0,
     .2,0,2,0,1,49*0,9,0,-9,0,-9,0,9,9*0,6,0,-6,0,3,0,-3,
     .0,6,0,3,0,-6,0,-3,0,4,0,2,0,2,0,1,0,-27,2*27,-27,27,
     .2*-27,27,-18,2*18,-18,-9,2*9,-9,-18,18,-9,9,18,-18,
     .9,-9,-18,-9,18,9,18,9,-18,-9,-12,12,-6,6,-6,6,-3,3,
     .-12,-6,12,6,-6,-3,6,3,-12,2*-6,-3,12,2*6,3,-8,2*-4,
     .-2,-4,2*-2,-1,18,2*-18,18,-18,2*18,-18,12,2*-12,12,
     .6,2*-6,6,12,-12,6,-6,-12,12,-6,6,2*9,4*-9,2*9,8,-8,
     .4,-4,4,-4,2,-2,2*6,2*-6,2*3,2*-3,2*6,2*3,2*-6,2*-3,
     .2*4,4*2,2*1,-6,0,6,0,6,0,-6,0,-4,0,4,0,-2,0,2,0,-3,
     .0,-3,0,3,0,3,9*0,-2,0,-2,0,-1,0,-1,49*0,-6,0,6,0,6,
     .0,-6,9*0,-4,0,4,0,-2,0,2,0,-3,0,-3,0,3,0,3,0,-2,0,
     .-2,0,-1,0,-1,0,18,2*-18,18,-18,2*18,-18,12,2*-12,12,
     .6,2*-6,6,9,-9,9,2*-9,9,-9,9,12,6,-12,-6,-12,-6,12,
     .2*6,-6,6,-6,3,-3,3,-3,8,4,-8,-4,4,2,-4,-2,6,3,6,3,
     .-6,-3,-6,-3,4,2,4,2*2,1,2,1,-12,2*12,-12,12,2*-12,
     .12,-8,2*8,-8,-4,2*4,-4,-6,6,-6,2*6,-6,6,3*-6,4*6,2*-6,
     .-4,4,-4,4,-2,2,-2,2,2*-4,2*4,2*-2,2*2,4*-3,4*3,4*-2,
     .4*-1,2,3*0,-2,3*0,1,3*0,1,75*0,2,3*0,-2,11*0,1,3*0,
     .1,19*0,-6,6,2*0,6,-6,2*0,-3,3,2*0,-3,3,10*0,-4,-2,
     .2*0,4,2,10*0,-2,-1,2*0,-2,-1,18*0,4,-4,2*0,-4,4,2*0,
     .2,-2,2*0,2,-2,10*0,2*2,2*0,2*-2,10*0,2*1,2*0,2*1,34*0,
     .2,3*0,-2,11*0,1,3*0,1,75*0,2,3*0,-2,3*0,1,3*0,1,19*0,
     .-6,6,2*0,6,-6,10*0,-3,3,2*0,-3,3,10*0,-4,-2,2*0,4,
     .2,2*0,-2,-1,2*0,-2,-1,18*0,4,-4,2*0,-4,4,10*0,2,-2,
     .2*0,2,-2,10*0,2*2,2*0,2*-2,2*0,2*1,2*0,2*1,2*0,-6,
     .0,6,0,6,0,-6,0,-3,0,3,0,-3,0,3,0,-4,0,-2,0,4,0,2,9*0,
     .-2,0,-1,0,-2,0,-1,49*0,-6,0,6,0,6,0,-6,9*0,-3,0,3,
     .0,-3,0,3,0,-4,0,-2,0,4,0,2,0,-2,0,-1,0,-2,0,-1,0,18,
     .2*-18,18,-18,2*18,-18,9,2*-9,2*9,2*-9,9,12,-12,6,-6,
     .-12,12,-6,6,12,6,-12,-6,-12,-6,12,2*6,-6,3,-3,6,-6,
     .3,-3,6,3,-6,-3,6,3,-6,-3,8,2*4,2,-8,2*-4,-2,4,2*2,
     .1,4,2*2,1,-12,2*12,-12,12,2*-12,12,-6,2*6,2*-6,2*6,
     .-6,-8,8,-4,4,8,-8,4,-4,2*-6,4*6,2*-6,-4,4,-2,2,-4,
     .4,-2,2,2*-3,2*3,2*-3,2*3,2*-4,2*-2,2*4,2*2,2*-2,2*-1,
     .2*-2,2*-1,4,0,-4,0,-4,0,4,0,2,0,-2,0,2,0,-2,0,2,0,
     .2,0,-2,0,-2,9*0,1,0,1,0,1,0,1,49*0,4,0,-4,0,-4,0,4,
     .9*0,2,0,-2,0,2,0,-2,0,2,0,2,0,-2,0,-2,0,1,0,1,0,1,
     .0,1,0,-12,2*12,-12,12,2*-12,12,-6,2*6,2*-6,2*6,2*-6,
     .6,-6,2*6,-6,6,-6,-8,-4,8,4,8,4,-8,-4,-3,3,-3,3,-3,
     .3,-3,3,-4,-2,4,2,-4,-2,4,2,-4,-2,-4,-2,4,2,4,2,-2,
     .-1,-2,-1,-2,-1,-2,-1,8,2*-8,8,-8,2*8,-8,4,2*-4,2*4,
     .2*-4,2*4,-4,4,2*-4,4,-4,3*4,4*-4,2*4,2,-2,2,-2,2,-2,
     .2,-2,2*2,2*-2,2*2,2*-2,4*2,4*-2,8*1/

c     set up dx, dy, dz:

      dx = xmax - xmin
      dy = ymax - ymin
      dz = zmax - zmin

c     pack the fxn value and derivative info for easier user later

c     since we first determine the coefficients on the unit cube,
c     which is (1/dx, 1/dy, 1/dz) times larger than the (dx,dy,dz)
c     box, we need to scale the derivatives appropriately.

      do 30 i = 0, 1
         do 20 j = 0, 1
            do 10 k = 0, 1
               finput(i, j, k, 1) = f(i, j, k)
               finput(i, j, k, 2) = fx(i, j, k) * dx
               finput(i, j, k, 3) = fy(i, j, k) * dy
               finput(i, j, k, 4) = fz(i, j, k) * dz
               finput(i, j, k, 5) = fxy(i, j, k) * dx * dy
               finput(i, j, k, 6) = fxz(i, j, k) * dx * dz
               finput(i, j, k, 7) = fyz(i, j, k) * dy * dz
               finput(i, j, k, 8) = fxyz(i, j, k) * dx * dy * dz
 10         continue
 20      continue
 30   continue

c     calculate the coefficients for a spline in the 
c     (0,0,0) x (1,1,1) cube. 
               
      do 160 l = 0, 3
         do 150 m = 0, 3
            do 140 n = 0, 3
               coeff(n, m, l) = 0.d0
               do 130 ideriv = 1,  8
                  do 120 i = 0, 1
                     do 110 j = 0, 1
                        do 100 k = 0, 1
                           coeff(n, m, l) = coeff(n, m, l)
     .                          + imat(k, j, i, ideriv, n, m, l)
     .                          * finput(i, j, k, ideriv)
 100                    continue
 110                 continue
 120              continue
 130           continue
 140        continue
 150     continue
 160  continue

c     at this point, spline coefficients are correct for the unit
c     cube, (0,0,0) x (1,1,1).  scale them so that they apply to
c     the box (0,0,0) x (dx,dy,dz)

c     f(x,y,z) = \sum_{l=0}^3 \sum_{m=0}^3 \sum_{n=0}^3 a_{lmn}
c                   ((x-xmin)/dx)^l ((y-ymin)/dx)^m ((z-zmin)/dz)^n
c              = \sum_{l=0}^3 \sum_{m=0}^3 \sum_{n=0}^3 
c                   a_{lmn} / dx^l / dy^m / dz^n 
c                   x (x-xmin)^l (y-ymin)^m (z-zmin)^n

c     precalculate dx^n, dy^n, dz^n:
      
      dxfac(0) = 1.d0
      dyfac(0) = 1.d0
      dzfac(0) = 1.d0
      do 165 i = 1, 3
         dxfac(i) = dxfac(i-1) * dx
         dyfac(i) = dyfac(i-1) * dy
         dzfac(i) = dzfac(i-1) * dz
 165  continue

      do 190 l = 0, 3
         do 180 m = 0, 3
            do 170 n = 0, 3
               coeff(n,m,l) = coeff(n,m,l) / dxfac(l) / dyfac(m)
     .              / dzfac(n)
 170        continue
 180     continue
 190  continue

c     convert spline coefficients to those correct for the 
c     (xmin,ymin,zmin) x (xmin+1,ymin+1,zmin+1) unit cube

c     f(x,y,z) = \sum_{l=0}^3 \sum_{m=0}^3 \sum_{n=0}^3 b_{lmn}
c                  (x-xmin)^l (y-ymin)^m (z-zmin)^n
c     using the binomial expansion,
c         = \sum_{l=0}^3 \sum_{m=0}^3 \sum_{n=0}^3 b_{lmn} \sum_{i=0}^l
c             (-1)^i xmin^i (l i) x^{l-i} \sum_{j=0}^m (-1)^j ymin^j
c             (m j) y^{m-j} \sum_{k=0}^n (-1)^k zmin^k (n k) z^{n-k}

c     zero the new coefficients, since they are calculated out of
c     order

      do 220 l = 0, 3
         do 210 m = 0, 3
            do 200 n = 0, 3
               coeff2(n,m,l) = 0.d0
 200        continue
 210     continue
 220  continue

c     precalculate xmin^n, ymin^n, zmin^n

      xmfac(0) = 1.d0
      ymfac(0) = 1.d0
      zmfac(0) = 1.d0
      do 225 i = 1, 3
         xmfac(i) = xmfac(i-1) * xmin
         ymfac(i) = ymfac(i-1) * ymin
         zmfac(i) = zmfac(i-1) * zmin
 225  continue

c     precalculate (-1)^n

      onefac(0) = 1.d0
      do 227 i = 1, 9
         onefac(i) = -onefac(i-1)
 227  continue

c     I think I could probably collapse this sextuple sum down to a
c     triple sum, but it's easier to let the computer do the work:

      do 280 l = 0, 3
         do 270 i = 0, l
            xfac = xmfac(i) * ichoos(l,i)
            do 260 m = 0, 3
               do 250 j = 0, m
                  yfac = ymfac(j) * ichoos(m,j)
                  do 240 n = 0, 3
                     do 230 k = 0, n
                        zfac = zmfac(k) * ichoos(n,k)
                        coeff2(n-k,m-j,l-i) = coeff2(n-k,m-j,l-i)
     .                       + onefac(i+j+k) * coeff(n,m,l) * xfac 
     .                       * yfac * zfac
 230                 continue
 240              continue
 250           continue
 260        continue
 270     continue
 280  continue

      return
      end
c
c-----------------------------------------------------------------------
c     qucof calculates the coefficients of a (one-dimensional) quintic
c     spline, given the value of the spline and its first two
c     derivatives at each of two knots.  Equations derived in SJS's
c     notebook on 12/20/01....sjs 12/21/01
c-----------------------------------------------------------------------
c
      subroutine qucof(x, y, dydx, d2ydx2, c)

      implicit real*8 (a-h,o-z)

c     x(1:2) contains the positions of the two spline knots, x1 and x2
c     y(1:2) contains the function values at the knots
c     dydx(1:2) contains the values of dy/dx at the knots
c     d2ydx2(1:2) contains the values of d2y/dx2 at the knots
c     c(1:6) contains the coefficients of the quintic spline
      dimension a(6), c(6),
     .     dydx(2), d2ydx2(2), x(2), y(2)

c     we are interested in the coefficients of a cubic spline,

c     y(x) = c(1) + c(2) * x + c(3) + x**2 + c(4) * x**3  + c(5) * x**4
c                 + c(6) * x**5

c     before calculating the c(i), we first calculate the coefficients
c     of a transformed quintic spline,

c     y(x) = f(t(x))
c     f(t) = a(1) + a(2) * x + a(3) * x**2 + a(4) * x**3 + a(5) * x**4
c                 + a(6) * x**5
c     t(x) = (x - x1) / (x2 - x1)

      delta = x(2) - x(1)
      a(1) = y(1)
      a(2) = delta * dydx(1)
      a(3) = 0.5d0 * delta**2 * d2ydx2(1)
      a(4) = -10.d0 * y(1) 
     .     - 6.d0 * delta * dydx(1) 
     .     - 3.d0 * 0.5d0 * delta**2 * d2ydx2(1)
     .     + 10.d0 * y(2) 
     .     - 4.d0 * delta * dydx(2)
     .     + 0.5d0 * delta**2 * d2ydx2(2)
      a(5) = 15.d0 * y(1) 
     .     + 8.d0 * delta * dydx(1)
     .     + 3.d0 * 0.5d0 * delta**2 * d2ydx2(1)
     .     - 15.d0 * y(2)
     .     + 7.d0 * delta * dydx(2)
     .     - 2.d0 * 0.5d0 * delta**2 * d2ydx2(2)
      a(6) = -6.d0 * y(1)
     .     - 3.d0 * delta * dydx(1)
     .     - 0.5d0 * delta**2 * d2ydx2(1)
     .     + 6.d0 * y(2)
     .     - 3.d0 * delta * dydx(2)
     .     + 0.5d0 * delta**2 * d2ydx2(2)

c     now calculate the coefficients of the y(x) spline using the
c     coefficients of the f(t) spline

      c(1) = a(1)
     .     - x(1) * a(2) / delta
     .     + x(1)**2 * a(3) / delta**2
     .     - x(1)**3 * a(4) / delta**3
     .     + x(1)**4 * a(5) / delta**4
     .     - x(1)**5 * a(6) / delta**5
      c(2) = a(2) / delta
     .     - 2.d0 * x(1) * a(3) / delta**2
     .     + 3.d0 * x(1)**2 * a(4) / delta**3
     .     - 4.d0 * x(1)**3 * a(5) / delta**4
     .     + 5.d0 * x(1)**4 * a(6) / delta**5
      c(3) = a(3) / delta**2
     .     - 3.d0 * x(1) * a(4) / delta**3
     .     + 6.d0 * x(1)**2 * a(5) / delta**4
     .     - 10.d0 * x(1)**3 * a(6) / delta**5
      c(4) = a(4) / delta**3
     .     - 4.d0 * x(1) * a(5) / delta**4
     .     + 10.d0 * x(1)**2 * a(6) / delta**5
      c(5) = a(5) / delta**4
     .     - 5.d0 * x(1) * a(6) / delta**5
      c(6) = a(6) / delta**5
      end
c
