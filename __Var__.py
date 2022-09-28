# Class for statistical variables, specifically serving for the quantification of source events,
# having two subclasses, VarSum and VarProd

import numpy as np

class Var:
    def __init__(self, value, error):
        self.val = value
        self.err = np.abs(error)
    def __add__(self, other):
        return Var(self.val + other, self.err)
    def __mul__(self, other):
        return Var(self.val * other, self.err * np.abs(other))
    def __pow__(self, other):
        return Var(self.val ** other, np.abs(other) * np.abs(self.val)**(other-1) * self.err)
    def __str__(self):
        stream = '\n'
        try:
            var_size = self.val.size
        except AttributeError:
            var_size = 1
        if var_size == 1:
            stream += ('{} \u00B1 {}\n').format(self.val, self.err)
        else:
            for i in np.arange(var_size):
                stream += ('{} \u00B1 {}\n').format(self.val[i], self.err[i])
        stream += '------------------------------------------------'
        return stream
    def sum_val(self):
        return Var(self.val.sum(),self.err.sum())

class VarSum(Var):
    # constructor taking var as a list of variables to be summed over and neq as positions for the minuses,
    # for example, A - B - C => neq = [1,2], position of A is 0
    def __init__(self, var, neg = []):
        self.val = 0
        sum_err = 0
        op = np.ones(len(var))
        for n in neg:
            op[n] = -1
        for i in np.arange(len(var)):
            self.val += var[i].val * op[i]
            sum_err += (var[i].err)**2
        self.err = np.sqrt(sum_err)

class VarProd(Var):
    # constructor taking var as a list of variables to be multiplied and inv as positions for the divisions,
    # for example, A / B / C => inv = [1,2], position of A is 0
    def __init__(self, var,  inv = []):
        self.val = 1
        sum_err = 0
        op = np.ones(len(var))
        for i in inv:
            op[i] = -1
        for i in np.arange(len(var)):
            self.val *= (var[i].val)**op[i]
            sum_err += (var[i].err/var[i].val)**2
        self.err = np.abs(self.val) * np.sqrt(sum_err)