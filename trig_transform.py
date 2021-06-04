

def trig_transform(self, x, u, debug=False):
    r"""
    Perform the trig_transform with simplification

    Supported functions:
        x = a * sin(t),
        x = a * cos(t),
        x = a / sin(t)
        x = a * tg(t)

    Examples
    ========

    >>> from sympy import Integral
    >>> from sympy import sqrt
    >>> from sympy.abc import x, t
    >>> p = Integral(sqrt(4 - x ** 2), (x, -2, -1))
    >>> p.trig_transform(x, 2 * sin(t))
    Integral(4*cos(t)**2, (t, -pi/2, -pi/6))

    Errors
    ======

    >>> from sympy import Integral
    >>> from sympy import sqrt
    >>> from sympy.abc import x, t
    >>> p = Integral(sqrt(4 - x ** 2), (x, -6, 7))
    >>> p.trig_transform(x, 2 * sin(t))

    Value error
    \not \exists t : 2 * sin(t) = -6

    Same value error for 2 / sin(t)

    Value error in cases:
        1. u contains more than 1 free variable
        2. u is not Symbol
        3. x is not Symbol

    See Also
    ========

    sympy.integrals.trigonometry.trigintegrate
    sympy.integrals.integrals.Integral.transform
    
    SymPy core methods used
    =======================

    sympy.solvers.solvers.solve
    sympy.core.basic.Basic.xreplace
    sympy.simplify.trigsimp
    """
    
    from sympy.solvers.solvers import solve, posify
    from sympy import cos, sin, tan, cot, sqrt

    d = Dummy('d')

    xfree = x.free_symbols.intersection(self.variables)
    if len(xfree) > 1:
        raise ValueError(
            'F(x) can only contain one of: %s' % self.variables)
    xvar = xfree.pop() if xfree else d

    if xvar not in self.variables:
        return self

    u = sympify(u)
    if isinstance(u, Expr):
        ufree = u.free_symbols
        if len(ufree) == 0:
            raise ValueError('f(u) cannot be a constant')
        if len(ufree) > 1:
            raise ValueError('More than 1 free symbol')
        uvar = ufree.pop()
    else:
        u, uvar = u

    if x.is_Symbol and u.is_Symbol:
        return self.xreplace({x: u})

    if not x.is_Symbol and not u.is_Symbol:
        raise ValueError('Either x or u must be a symbol')

    if uvar == xvar:
        return self.transform(x, (u.subs(uvar, d), d)).xreplace({d: uvar})

    if uvar in self.limits:
        raise ValueError('u must contain the same variable as in x or a variable that is not already an integration variable')

    A = u.args[0] if len(u.args) > 1 else 1

    f = self.args[0]
    newF = f.subs(x, u)

    if debug:
        print("New f:", newF)

    dfAndlimits = self.args[1]
    newDf = self.args[1][0].subs(x, u)

    if debug:
        print("New df:", newDf)

    if len(dfAndlimits) == 1:  # indefinite integral
        newI = self.func(newF * diff(newDf))

        return newI.trigsimp().xreplace({
            Abs(cos(uvar)): cos(uvar),
            Abs(sin(uvar)): sin(uvar),
            Abs(tan(uvar)): tan(uvar),
            Abs(cot(uvar)): cot(uvar)
        }).trigsimp()

    # definite integral
    def __calc_limits(a, b):
        funcs = solve(u - x, uvar)
        F = funcs[-1]

        # exception if wrong limits for sin
        if u.args[-1] == sin(uvar) and (abs(a / A) > 1 or abs(b / A) > 1):
            raise ValueError("Wrong limits")


        return F.subs(xvar, a), F.subs(xvar, b)

    oldLimits = (self.limits[0][1], self.limits[0][2])

    b, a = __calc_limits(oldLimits[0], oldLimits[1])
    
    c = 1
    
    if a < b:
        a, b = b, a
        c = -1
    
    newI = self.func(c * newF * diff(newDf), (uvar, b, a))

    if debug:
        print("New limits:", b, a)

    return newI.trigsimp().xreplace({
        Abs(cos(uvar)): cos(uvar),
        Abs(sin(uvar)): sin(uvar),
        Abs(tan(uvar)): tan(uvar),
        Abs(cot(uvar)): cot(uvar)
    }).trigsimp()
