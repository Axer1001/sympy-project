from sympy import Symbol
from sympy import sqrt, sin, cos, tan, atan
from sympy import Integral
from sympy import Abs

x = Symbol('x')
t = Symbol('t', real=True)

i1 = Integral(sqrt(4-x**2), (x, -2 , -1))
i1 = i1.trig_transform(x, 2 * sin(t))

i2 = Integral(sqrt(x**2-4), (x, 3, 4))
i2 = i2.trig_transform(x, 2/sin(t))

i3 = Integral(sqrt(x**2 + 4), (x, -3, 1))
i3 = i3.trig_transform(x, 2 * tan(t))

# sin(x) = 5 exception!
#i1 = i1.xreplace({Abs(cos(t)) : -cos(t)})
print(i1.doit(), i2.doit(), i3.doit(), sep='\n')



def trig_transform(self, x, u, debug=False, t=0):
    
    def empty_handler(signal, f):
        return

    if t:
        import signal

        def signal_handler(signal, f):
            raise Exception("Time Out!")
    
        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(t)

    from sympy.solvers.solvers import solve, posify
    from sympy import cos, sin, tan, cot, sqrt

    d = Dummy('d')

    xfree = x.free_symbols.intersection(self.variables)
    if len(xfree) > 1:
        signal.signal(signal.SIGALRM, empty_handler)
        raise ValueError(
            'F(x) can only contain one of: %s' % self.variables)
    xvar = xfree.pop() if xfree else d

    if xvar not in self.variables:
        return self

    u = sympify(u)
    if isinstance(u, Expr):
        ufree = u.free_symbols
        if len(ufree) == 0:
            signal.signal(signal.SIGALRM, empty_handler)
            raise ValueError(filldedent('''
                f(u) cannot be a constant'''))
        if len(ufree) > 1:
            signal.signal(signal.SIGALRM, empty_handler)
            raise ValueError(filldedent('''
                When f(u) has more than one free symbol, the one replacing x
                must be identified: pass f(u) as (f(u), u)'''))
        uvar = ufree.pop()
    else:
        u, uvar = u
        if uvar not in u.free_symbols:
            signal.signal(signal.SIGALRM, empty_handler)
            raise ValueError(filldedent('''
                Expecting a tuple (expr, symbol) where symbol identified
                a free symbol in expr, but symbol is not in expr's free
                symbols.'''))

    if x.is_Symbol and u.is_Symbol:
        signal.signal(signal.SIGALRM, empty_handler)
        return self.xreplace({x: u})

    if not x.is_Symbol and not u.is_Symbol:
        signal.signal(signal.SIGALRM, empty_handler)
        raise ValueError('either x or u must be a symbol')

    if uvar == xvar:
        signal.signal(signal.SIGALRM, empty_handler)
        return self.transform(x, (u.subs(uvar, d), d)).xreplace({d: uvar})

    if uvar in self.limits:
        signal.signal(signal.SIGALRM, empty_handler)
        raise ValueError(filldedent('''
            u must contain the same variable as in x
            or a variable that is not already an integration variable'''))

    A = u.args[0] if len(u.args) > 1 else 1

    f = self.args[0]
    newF = f.subs(x, u)

    if debug:
        print("New f:", newF)

    dfAndlimits = self.args[1]
    newDf = self.args[1][0].subs(x, u)

    if debug:
        print("New df:", newDf)

    if len(dfAndlimits) == 1: # indefinite integral
        newI = self.func(newF * diff(newDf))
        
        signal.signal(signal.SIGALRM, empty_handler)
        return newI.trigsimp().xreplace({
            Abs(cos(uvar)) : cos(uvar),
            Abs(sin(uvar)) : sin(uvar),
            Abs(tan(uvar)) : tan(uvar),
            Abs(cot(uvar)) : cot(uvar)
        }).trigsimp()

    # definite integral
    def __calc_limits(a, b):
        funcs = solve(u - x, uvar)
        F = funcs[-1]

        # exception if wrong limits for sin
        if u.args[-1] == sin(uvar) and (abs(a / A) > 1 or abs(b / A) > 1):
            signal.signal(signal.SIGALRM, empty_handler)
            raise(ValueError("Wrong limits"))
        
        signal.signal(signal.SIGALRM, empty_handler)
        return F.subs(xvar, a), F.subs(xvar, b)

    oldLimits = (self.limits[0][1], self.limits[0][2])

    a, b = __calc_limits(oldLimits[0], oldLimits[1])
    newI = self.func(newF * diff(newDf), (uvar, a, b))
    
    if debug:
        print("New limits:", a, b)

    signal.signal(signal.SIGALRM, empty_handler)
    return newI.trigsimp().xreplace({
        Abs(cos(uvar)) : cos(uvar),
        Abs(sin(uvar)) : sin(uvar),
        Abs(tan(uvar)) : tan(uvar),
        Abs(cot(uvar)) : cot(uvar)
    }).trigsimp()
