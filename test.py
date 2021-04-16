from sympy.abc import a, x, u
from sympy import Integral, cos, sin, tan
from sympy import trigsimp
i = Integral((4 - x ** 2)**(1/2), x)
#i = Integral(1/(2 * sin(x) + 3 * cos(x) - 1), x)
i = i.trig_transform(x, 2 * sin(a))
print(i.trigsimp())
