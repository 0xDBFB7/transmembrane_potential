
"""

MOVED TO TRANSMEMBRANE_MATH.IPYNB

THIS IS NOW FIXED - switched A and B coefficients appropriately

PSOPT_Manual_R5.pdf

For convenience, let's convert the Kotnik equation to standard DE form.

per https://lpsa.swarthmore.edu/Representations/SysRepTransformations/TF2SDE.html

originally (mistakenly) used A6c (1998). Switched to eq 8, (multiplied by R, eq 10).
also missed switching the top and bottom: the top goes with the input terms (right) and the bottom with the output.

H(s)= (R*X(s)) / U(s) = (R a3 s^2 + R a2 s + R a3) / (b3 s^2 + b2 s + b1)

"Solution: Separate the equation so that the output terms, X(s), are on the
left and the input terms, Fa(s), are on the right.  Make sure there are only positive powers of s."

(b3 s^2 + b2 s + b1) X(s) = (R a3 s^2 + R a2 s + R a3) U(s)

"Now take the inverse Laplace Transform (so multiplications by "s" in the Laplace domain are replaced by derivatives in time)."

b3 x'' + b2 x' + b1 x = R a3 u'' + R a2 u' + R a3 u

where u is the input function,

now, gekko wants it in first order, so we have to substitute (Convert to a System of DEs)
http://www.math.utah.edu/~gustafso/2250systems-de.pdf
http://www.sharetechnote.com/html/DE_HigherOrderDEtoFirstOrderDE.html
https://math.berkeley.edu/~zworski/128/psol12.pdf

https://math.stackexchange.com/questions/1120984/

(reworked because I found the notation confusing (u3 should be the third derivative!))

u0 = u
u1 = u'
u2 = u''
u3 = u'''

      x0 = x
x0' = x1 = x'
x1' = x2 = x'' =

    b3 x'' + b2 x' + b1 x = R a3 u'' + R a2 u' + R a3 u
    x''  = (R a3 u'' + R a2 u' + R a3 u - b2 x' - b1 x)/b3
    x2  = (R*a3*u2 + R*a2*u1 + R*a3*u0 - b2*x1 - b1*x0)/b3

docs: "In all simulation modes (IMODE=1,4,7), the number of equations must equal the number of variables."







"""
