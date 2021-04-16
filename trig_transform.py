def trig_transform(self, x, u):
    from sympy.core.sympify import sympify
    from sympy.core.symbol import Dummy, Symbol
    from sympy.core.expr import Expr
    from sympy import filldedent
    from sympy.series import limit
    from sympy.solvers.solvers import solve, posify
    d = Dummy('d')

    freeSymbols = x.free_symbols.intersection(self.variables)
    if len(freeSymbols) > 1:
        raise ValueError('Multiple variables')
    variable = freeSymbols.pop() if freeSymbols else d

    if variable not in self.variables:
        return self

    u = sympify(u) # transform to sympy type
    if isinstance(u, Expr):
        ufree = u.free_symbols
        if len(ufree) == 0:
            raise ValueError('Constant')
        if len(ufree) > 1:
            raise ValueError('Multiple variables')
        uvar = ufree.pop()
    else:
        u, uvar = u
        if uvar not in u.free_symbols:
            raise ValueError('Wrong parameters')
    if x.is_Symbol and u.is_Symbol:
        return self.xreplace({x: u})

    if not x.is_Symbol and not u.is_Symbol:
        raise ValueError('either x or u must be a symbol')

    if uvar == variable:
        return self.transform(x, (u.subs(uvar, d), d)).xreplace({d: uvar})

    if uvar in self.limits:
        raise ValueError('Variables error')

    if not x.is_Symbol:
        F = [x.subs(variable, d)]
        soln = solve(u - x, variable, check=False)
        if not soln:
            raise ValueError('no solution for solve(F(x) - f(u), x)')
        f = [fi.subs(uvar, d) for fi in soln]
    else:
        f = [u.subs(uvar, d)]
        pdiff, reps = posify(u - x)
        puvar = uvar.subs([(v, k) for k, v in reps.items()])
        soln = [s.subs(reps) for s in solve(pdiff, puvar)]
        if not soln:
            raise ValueError('no solution for solve(F(x) - f(u), u)')
        F = [fi.subs(variable, d) for fi in soln]

    newfuncs = {(self.function.subs(variable, fi) * fi.diff(d)
                 ).subs(d, uvar) for fi in f}
    newfunc = newfuncs.pop()

    def _calc_limits(limit1, limit2):
        # need to implement
        return limit1, limit2

    newlimits = _calc_limits(0, 0) # template
    return self.func(newfunc, *newlimits)

