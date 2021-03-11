from CoolProp.CoolProp import PropsSI, PhaseSI

class Medium():
    """ 
    Class Medium slouzi k uchovani informaci o mediich vstupujicich do vymeniku.
    Ve tride jsou dale vypocitany fyzikalni vlasnosti nutne k dalsim vypoctum.
     """
    def __init__(self, part : dict, power : float, hotter : bool):
        """ 
        Parameters
        ----------

        part : dict
            Obsahuje veskere informace o mediu. Teplota, tlak, slozeni, zanaseni a prutok.
        power : float
            Informace o vykonu ve Wattech.
        hotter : bool
            Informace jestli se bude medium chladit nebot oteplovat.

        ----------
         """
        self.t1 = part['T1']
        self.p = part['P']
        self.medium = part['Medium']
        self.zanaseni = part['Rf']
        self.prutok = part['M']
        self.mol = self.molar()
        self.t2 = self.temp_out(power, hotter)
        self.t_str = (self.t2 + self.t1) / 2
        self.props(hotter)

    def props(self, hotter:bool):
        """ 
        Metoda slouzi k dopocitani dalsich fyzikalnich vlastnosti medii
         """
        self.rho = self.propDMASS()
        self.cp = self.propCPMASS()
        self.eta = self.propViscosity()
        self.lamb = self.propConductivity()
        self.phase = self.propPhase()
        if self.phase == 'twophase': raise Exception('V {} mediu dochazi ke zmene skupenstvi.'.format('teplejsim' if hotter else 'chladnejsim'))
    def temp_out(self, power : float, hotter : bool) -> float:
        """
        Metoda slouzi k vypocitani vystupni teploty. Pro vypocet je pouzita metoda puleni intervalu.
        ----------

        power : float
            Informace o vykonu ve Wattech.
        hotter : bool
            Informace jestli se bude medium chladit nebot oteplovat.

        ----------        
        """
        h1 = self.propH(self.t1)
        if hotter:
            h2 = h1 - power / self.prutok
            t_1 = 274
            t_2 = self.t1    
        else:
            h2 = h1 + power / self.prutok
            t_1 = self.t1
            t_2 = 2000
        while 1:
            temp = (t_1 + t_2)/2
            h2_teo = self.propH(temp)
            diff = h2 - h2_teo
            if abs(diff) < 10:
                break
            elif diff > 0:
                t_1 = temp
            else:
                t_2 = temp
        return temp
    def molar(self):
        """ 
        Metoda slouzi ke spocitani molarnich hmotnosti pro jednotliva media
         """
        M = []
        for i in self.medium:
            M.append(PropsSI('MOLARMASS', 'T', self.t1, 'P', self.p, i[0]))
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
    def propH(self, temperature : float):
        """ 
        Metoda slouzi k vypoctu rozdilu entaplii na zaklade teploty pro jednotliva media
        ----------

        temperature : float
            Teplota pro kterou ma byt entalpie vypocitana.

        ----------        
        """
        a = 0
        b = 0
        for i in range(len(self.medium)):
            if self.medium[i][0] == 'H2O':
                enthalpy  = PropsSI('HMASS', 'T', temperature, 'P', self.p*self.medium[i][1], 'H2O')
            else:    
                enthalpy  = PropsSI('HMASS', 'T', temperature, 'P', self.p, self.medium[i][0])
            a += self.medium[i][1] * self.mol[i] * enthalpy 
            b += self.medium[i][1] * self.mol[i]
        return a/b
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
    '''
    medium_vzor1 = {'M': 11.186/3600, 'T1': 900+273.15, 'T2': 65.81+ 273.15, 'P': 100000, 'Rf': 10.0, 'Medium': [['O2', 0.07], ['CO2', 0.12], ['Ar', 0.006], ['H2O', 0.1], ['N2', 0.7]]}
    medium_vzor2 = {'M': 0.02894, 'T1': 288.15, 'T2': 313.15, 'P': 500000.0, 'Rf': 10.0, 'Medium': [['Water', 1]]}
    '''
    medium_vzor1 = {'M': 150/3600, 'T1': 890+273.15, 'P': 100000, 'Rf': 10.0, 'Medium': [['O2', 0.07], ['CO2', 0.12], ['Ar', 0.006], ['H2O', 0.1], ['N2', 0.7]]}
    medium_vzor2 = {'M': 1.9279, 'T1': 298.15, 'P': 400000.0, 'Rf': 10.0, 'Medium': [['Water', 1]]}
    medium1 = Medium(medium_vzor1, 40290, True)
    medium2 = Medium(medium_vzor2, 40290, False)
    print(medium2.phase)
    print(medium2.t1,medium2.t2)
    print(medium1.phase)
    print(medium1.t1,medium1.t2)
    
    print('SUCCES')