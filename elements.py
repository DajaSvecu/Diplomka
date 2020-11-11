from CoolProp.CoolProp import PropsSI
class Element(object):
    def __init__(self, x, M):
        self.x = x
        self.M = M
        self.slozky = ['O2', 'CO2', 'N2', 'Ar', 'H2O']

    def w(self):
        # VYPOCET HMOTNOSTNI ZASTOUPENI PRVKU
        a = []
        for i in range(5):
            a.append(self.x[i]*self.M[i])
        b = []
        for i in range(5,0,-1):
            b.append(a[-i]/(a[-i] + a[-i+1] + a[-i+2] + a[-i+3] + a[-i+4]))
        return b
    def enthalpy(self, t, p):
        h = [0,0,0,0,0]
        for i in range(4):
            h[i] = PropsSI('H', 'T', t, 'P', p, self.slozky[i]) - PropsSI('H', 'T', 273.15, 'P', p, self.slozky[i])
        h[-1] = PropsSI('H', 'T', t, 'P', p * self.x[-1], 'Water') - PropsSI('H', 'T', 273.16, 'P', p * self.x[-1], 'Water')
        zlom = self.w()
        h_v = 0
        for i in range(5):
            h_v += h[i]*zlom[i]
        return h_v
    def enthalpy_w(self, t, p):
        return PropsSI('H', 'T', t, 'P', p, 'Water')