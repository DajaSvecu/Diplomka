from CoolProp.CoolProp import PropsSI, PhaseSI
def props(part:dict) -> None:
    t_str = (part['T1'] + part['T2'])/2
    part['deltaH'] = propH(part)
    for name in ['Dmass', 'Cpmass', 'viscosity', 'conductivity', 'Prandtl']:
        part[name] = prop(name, part, t_str)
    part['Q'] = part['M'] * part['deltaH']
def propH(medium:dict) -> float:
        """ 
        Metoda slouzi ke spocitani rozdilu entalpie na zaklade teplot
         """
        return abs(PropsSI('HMASS', 'T', medium['T1'], 'P', medium['P'], medium['Medium']) - PropsSI('HMASS', 'T', medium['T2'], 'P', medium['P'], medium['Medium']))
def prop(name : str, medium : dict, t : float) -> float:
        """ 
        Metoda slouzi ke spocitani dane fyzikalni vlastnosti
        ----------
        PARAMETERS
        ----------
        name : str
            Nazev fyzikalni vlastnosti
         """
        return PropsSI(name, 'T', t, 'P', medium['P'], medium['Medium'])

class Medium():
    """ 
    Class Medium slouzi k uchovani informaci o mediich vstupujicich do vymeniku.
    Ve tride jsou dale vypocitany fyzikalni vlasnosti nutne k dalsim vypoctum.
     """
    def __init__(self, part : dict):
        """ 
        Parameters
        ----------

        part : dict
            Obsahuje veskere informace o mediu. Teplota, tlak, slozeni, zanaseni a prutok.
        
        ----------
         """
        self.t1 = part['T1']
        self.t2 = part['T2']
        self.p = part['P']
        self.medium = part['Medium']
        self.zanaseni = part['Rf']
        self.prutok = part['M']
        self.t_str = (self.t2 + self.t1) / 2
        self.props()
    def props(self):
        """ 
        Metoda slouzi k dopocitani dalsich fyzikalnich vlastnosti medii
         """
        self.mol = self.molar()
        self.rho = self.propDMASS()
        self.cp = self.propCPMASS()
        self.eta = self.propViscosity()
        self.lamb = self.propConductivity()
        self.deltaH = self.propH()
        self.phase = self.propPhase()
        self.q = self.deltaH * self.prutok
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
            a += self.medium[i][1] * self.mol[i] * PropsSI('CPMASS', 'T', self.t_str, 'P', self.p, self.medium[i][0])
            b += self.medium[i][1] * self.mol[i]
        return a/b
    def propViscosity(self):
        """ 
        Metoda slouzi k vypoctu viskozity media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
            a += self.medium[i][1] * self.mol[i]**0.5 * PropsSI('VISCOSITY', 'T', self.t_str, 'P', self.p, self.medium[i][0]) 
            b += self.medium[i][1] * self.mol[i] **0.5
        return a/b       
    def propConductivity(self):
        """ 
        Metoda slouzi k vypoctu tepelne vodivosti media
         """
        a = 0
        b = 0
        for i in range(len(self.medium)):
                a += self.medium[i][1] / self.mol[i]**0.5 * PropsSI('CONDUCTIVITY', 'T', self.t_str, 'P', self.p, self.medium[i][0]) 
                b += self.medium[i][1] / self.mol[i] **0.5
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
            a += self.medium[i][1] * self.mol[i] * difference
            b += self.medium[i][1] * self.mol[i]
        return abs(a/b)
    def propPhase(self):
        phases = []
        for i in self.medium:
            phase1 = PhaseSI('T', self.t1, 'P', self.p*i[1], i[0])
            if phase1 == 'supercritical_gas':
                phase1 = 'gas'
            elif phase1 == 'supercritical_liquid':
                phase1 = 'liquid'
            
            phase2 = PhaseSI('T', self.t2, 'P', self.p*i[1], i[0])
            if phase2 == 'supercritical_gas':
                phase2 = 'gas'
            elif phase2 == 'supercritical_liquid':
                phase2 = 'liquid'
            
            if phase1 == phase2:
                phases.append(phase1)
            else:
                return 'twophase'

        for phase in phases:
            if phase != phases[0]:
                return 'twophase'

        return phases[0]

if __name__ == '__main__':
    medium1 = {'M': 0.833, 'T1': 623.15, 'T2': 423.15, 'P': 100000.0, 'Rf': 10.0, 'Medium': 'Water'}
    multimedium_vzor = {'M': 2, 'T1': 298.15, 'T2': 348.15, 'P': 400000.0, 'Rf': 10.0, 'Medium': [['Water', 1]]}
    props(medium1)
    medium2 = Medium(multimedium_vzor)
    print(medium2.phase)