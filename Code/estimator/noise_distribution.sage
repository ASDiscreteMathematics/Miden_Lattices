class NoiseDistributionWithOverflow:
    def __init__( self, bound, cutoff_point, dictionary ):
        self.bound = bound
        self.cutoff_point = cutoff_point
        self.dictionary = dictionary

    def __add__( self, other ):
        RR = self.dictionary[0].parent()
        assert(self.cutoff_point == other.cutoff_point)
        bound = self.bound + other.bound

        if len(self.dictionary) == 1:
            return other

        left_max = max(self.dictionary)
        right_max = max(other.dictionary)
        left_list = [self.dictionary.get(n, RR(0)) for n in range(-left_max, left_max+1)]
        right_list = [other.dictionary.get(n, RR(0)) for n in range(-right_max, right_max+1)]

        conv = convolution(left_list, right_list)
        lenconv = len(conv)

        dictionary = dict()
        for n in range(lenconv):
            index = n - (lenconv-1)/2
            if index < -cutoff_point or index > cutoff_point:
                continue
            dictionary[index] = conv[n]

        return NoiseDistributionWithOverflow(bound, self.cutoff_point, dictionary)

    def __mul__( self, other ):
        RR = self.dictionary[0].parent()
        if str(type(other)) == "<class 'sage.rings.integer.Integer'>":
            acc = NoiseDistributionWithOverflow(self.bound, self.cutoff_point, {0:RR(1)})
            for b in bin(other)[2:]:
                acc = acc + acc
                if b == '1':
                    acc = acc + self
            return acc

        elif str(type(other)) == "<class '__main__.NoiseDistributionWithOverflow'>":
            assert(self.cutoff_point == other.cutoff_point)
            bound = self.bound * other.bound
            dictionary = dict()
            for l in self.dictionary:
                for r in other.dictionary:
                    if l * r < -self.cutoff_point or l * r > self.cutoff_point:
                        continue
                    else:
                        dictionary[r*l] = dictionary.get(r*l, RR(0)) + self.dictionary[l] * other.dictionary[r]
            return NoiseDistributionWithOverflow(bound, self.cutoff_point, dictionary)

    def var( self ):
        RR = self.dictionary[0].parent()
        acc = RR(0)
        for n in self.dictionary:
            acc += RR(n^2) * self.dictionary[n]
        E = self.E()
        var_lb = acc - E
        var_ub = acc + self.overflow() * self.bound^2 - E
        return var_lb, var_ub

    def std( self ):
        var_lb, var_ub = self.var()
        return sqrt(var_lb), sqrt(var_ub)

    def E( self ):
        # ignores the overflow term
        RR = self.dictionary[0].parent()
        acc = RR(0)
        for n in self.dictionary:
            acc += RR(n) * self.dictionary[n]
        return acc

    def overflow( self ):
        RR = self.dictionary[0].parent()
        acc = RR(0)
        for d in self.dictionary:
            acc += self.dictionary[d]
        return RR(1) - acc

    def statistical_distance( self, other ):
        # AKA total variational distance (TVD)
        assert(self.cutoff_point == other.cutoff_point), f"cutoff points: {self.cutoff_point} versus {other.cutoff_point}"
        RR = self.dictionary[0].parent()
        acc = RR(0)
        for n in range(-self.cutoff_point, self.cutoff_point+1):
            acc += abs(self.dictionary.get(n, RR(0)) - other.dictionary.get(n, RR(0)))
        return acc / RR(2), (acc + self.overflow() + other.overflow()) / RR(2)

    @staticmethod
    def CenteredBinomial(cutoff_point, precision, ncrumbs):
        RR = Reals(precision)
        binomial_dictionary = dict()
        for i in range(ncrumbs+1):
            binomial_dictionary[i] = binomial_dictionary.get(i, RR(0)) + binomial(ncrumbs,i) * RR(0.5)^ncrumbs
        centered_binomial_dictionary = dict()
        for bdl in binomial_dictionary:
            for bdr in binomial_dictionary:
                centered_binomial_dictionary[bdl-bdr] = centered_binomial_dictionary.get(bdl-bdr, RR(0)) + binomial_dictionary[bdl] * binomial_dictionary[bdr]
        for cbd in centered_binomial_dictionary:
            if cbd < -cutoff_point or cbd > cutoff_point:
                centered_binomial_dictionary.remove(cbd)
        return NoiseDistributionWithOverflow(ncrumbs, cutoff_point, centered_binomial_dictionary)

    @staticmethod
    def DiscreteGaussian(cutoff_point, precision, sigma):
        RR = Reals(precision)
        bound = infinity
        dictionary = dict()
        for n in range(-cutoff_point, cutoff_point+1):
            P = RR(1)/(RR(sigma) * sqrt(RR(2) * RR(pi))) * exp(-RR(0.5) * (n / RR(sigma))^2)
            dictionary[n] = P
        return NoiseDistributionWithOverflow(bound, cutoff_point, dictionary)
