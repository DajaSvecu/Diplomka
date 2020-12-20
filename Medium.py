from CoolProp.CoolProp import PropsSI
def props(medium:dict) -> None:
    t_str = (medium['T1'] + medium['T2'])/2
    medium['deltaH'] = propH(medium)
    for name in ['Dmass', 'Cpmass', 'viscosity', 'conductivity', 'Prandtl']:
        medium1[name] = prop(name, medium, t_str)
def propH(medium:dict) -> float:
        """ 
        Metoda slouzi ke spocitani rozdilu entalpie na zaklade teplot
         """
        return abs(PropsSI('HMASS', 'T', medium['T1'], 'P', medium['P'], medium['Medium']) - PropsSI('HMASS', 'T', medium['T2'], 'P', medium['P'], medium['Medium']))

def prop(name:str, medium : dict, t : float) -> float:
        """ 
        Metoda slouzi ke spocitani dane fyzikalni vlastnosti
        ----------
        PARAMETERS
        ----------
        name : str
            Nazev fyzikalni vlastnosti
         """
        return PropsSI(name, 'T', t, 'P', medium['P'], medium['Medium'])

if __name__ == '__main__':
    medium1 = {'M': 1, 'T1': 500, 'T2': 400, 'P': 100000.0, 'Rf': 10.0, 'Medium': 'Water'}
    props(medium1)




class MultiMedium():
    """ 
    Class MultiMedium slouzi k uchovani informaci o mediich vstupujicich do vymeniku.
    Ve tride jsou dale vypocitany fyzikalni vlasnosti nutne k dalsim vypoctum.
     """
    def __init__(self, t1, t2, p, medium):
        """ 
        Parameters
        ----------
        t1 : int/float 
            Vstupni teplota [K]
        t2 : int/float
            Vystupni teplota [K]
        p : int/float
            Vstupni tlak [Pa]
        medium : list
            Smes medii a jejich zastoupeni
         """
        self.t1 = t1
        self.t2 = t2
        self.p = p
        self.medium = medium
        self.t_str = (self.t2 + self.t1) / 2
    def props(self):
        """ 
        Metoda slouzi k dopocitani dalsich fyzikalnich vlastnosti medii
         """
        self.M = self.molar()
        self.rho = self.propDMASS()
        self.cp = self.propCPMASS()
        self.eta = self.propViscosity()
        self.lamb = self.propConductivity()
        self.deltaH = self.propH()

    def molar(self):
        """ 
        Metoda slouzi ke spocitani molarnich hmotnosti pro jednotliva media
         """
        M = []
        for i in self.medium:
            M.append(PropsSI('MOLARMASS', 'T', self.t_str, 'P', self.p, i[0]))
        return M
    def propDMASS(self):
        """ 
        Metoda slouzi k vypoctu hustoty media
         """
        a = 0
        for i in self.medium:
            a += i[1] * PropsSI('DMASS', 'T', self.t_str, 'P', self.p, i[0])
        return a

    def propCPMASS(self):
        """ 
        Metoda slouzi k vypoctu tepelne kapacity media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
            a += self.medium[i][1] * self.M[i] * PropsSI('CPMASS', 'T', self.t_str, 'P', self.p, self.medium[i][0])
            b += self.medium[i][1] * self.M[i]
        return a/b

    def propViscosity(self):
        """ 
        Metoda slouzi k vypoctu viskozity media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
            a += self.medium[i][1] * self.M[i]**0.5 * PropsSI('VISCOSITY', 'T', self.t_str, 'P', self.p, self.medium[i][0]) 
            b += self.medium[i][1] * self.M[i] **0.5
        return a/b       

    def propConductivity(self):
        """ 
        Metoda slouzi k vypoctu tepelne vodivosti media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
                a += self.medium[i][1] / self.M[i]**0.5 * PropsSI('CONDUCTIVITY', 'T', self.t_str, 'P', self.p, self.medium[i][0]) 
                b += self.medium[i][1] / self.M[i] **0.5
        return a/b 
    def propH(self):
        """ 
        Metoda slouzi k vypoctu rozdilu entaplii na zaklade teploty pro jednotliva media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
            if self.medium[i][0] == 'H2O':
                difference = PropsSI('HMASS', 'T', self.t1, 'P', self.p*self.medium[i][1], 'H2O') - PropsSI('HMASS', 'T', self.t2, 'P', self.p*self.medium[i][1], 'H2O')
            else:    
                difference = PropsSI('HMASS', 'T', self.t1, 'P', self.p, self.medium[i][0]) - PropsSI('HMASS', 'T', self.t2, 'P', self.p, self.medium[i][0])
            a += self.medium[i][1] * self.M[i] * difference
            b += self.medium[i][1] * self.M[i]
        return abs(a/b)  