from Medium import Medium
from heat_exchanger_inputs import Sizes
import math

class Calculate():
    def __init__(self, tube:dict, shell:dict, rest:dict):
        self.tube = Medium(tube) # TODO error handling two_phase state
        self.shell = Medium(shell) # TODO error handling two_phase state
        self.rest = rest
        self.sizes = Sizes()
        if self.tube.q > 0.99 * self.shell.q and self.tube.q < 1.01 * self.shell.q:
            self.rest['Q'] = 1.08 * (self.tube.q + self.shell.q) / 2
        #else: TODO error handling wrong power
    
    def calculate_all(self) -> None:
        heat_exchanger = {}
        vysledky = []
        pocet = 0
        pocet_exc = 0
        pocet_ass = 0
        for shell in self.sizes.list_of_shells:
            if shell['D2'] > self.rest['MaxSirka']: continue # podminka pro maximalni sirku
            heat_exchanger['shell'] = shell
            for tube in self.sizes.list_of_tubes:
                heat_exchanger['d_in'] = tube['d_in']
                for tl in tube['wall']:
                    heat_exchanger['tl'] = tl
                    for t_p in {0.2, 0.4, 0.6, 0.8, 1}:
                        heat_exchanger['t_p'] = round(t_p * (heat_exchanger['shell']['D2'] - 2 * heat_exchanger['shell']['wall']), 5)
                        # TODO pridelat ruzne vysky prepazky
                        for t_t in {1.25, 1.33, 1.5}:
                            heat_exchanger['t_t'] = round(t_t * (heat_exchanger['d_in'] + 2 * heat_exchanger['tl']), 5)
                            try:
                                heat_exchanger_result = self.calculate_heat_exchanger(heat_exchanger)
                            except AssertionError:
                                pocet_ass += 1
                            except Exception:
                                pocet_exc += 1
                            else:
                                vysledky.append(heat_exchanger_result)
                                pocet += 1
        print('TLAKOVA ZTRATA')
        vysledky.sort(key=lambda a: a[1]['tlak_ztraty'])
        for i in range(10):
            print(vysledky[i][1])

        print('VAHA')
        vysledky.sort(key=lambda a: a[1]['hmotnost'])
        for i in range(10):
            print(vysledky[i][1])
        
        print('KOMPAKTNOST')
        vysledky.sort(reverse=True, key=lambda a: a[1]['kompaktnost'])
        for i in range(10):
            print(vysledky[i][1])


        print(pocet_ass)
        print(pocet_exc)
        print(pocet)

    def calculate_heat_exchanger(self, heat_exchanger:dict) -> dict:
        D1 = heat_exchanger['shell']['D2'] - 2 * heat_exchanger['shell']['wall']
        DS = D1 - 0.016
        d_out = heat_exchanger['d_in'] + 2 * heat_exchanger['tl']
        h_p = 0.75 * D1 # ULOZIT PARAMETER
        s_p = self.baffle_thickness(D1, heat_exchanger['t_p']) # ULOZIT PARAMETER


        if self.rest['Uhel'] == 30:
            t_t1 = heat_exchanger['t_t'] * 1 
            t_t2 = heat_exchanger['t_t'] * 0.866
        elif self.rest['Uhel'] == 45:
            t_t1 = heat_exchanger['t_t'] * 1.414
            t_t2 = heat_exchanger['t_t'] * 0.707
        elif self.rest['Uhel'] == 60:
            t_t1 = heat_exchanger['t_t'] * 1.732
            t_t2 = heat_exchanger['t_t'] * 0.5
        elif self.rest['Uhel'] == 90:
            t_t1 = heat_exchanger['t_t'] * 1
            t_t2 = heat_exchanger['t_t'] * 1
        else:
            raise Exception('')

        # -----------------------TRUBKOVY PROSTOR-----------------------
        b1 = 0.866 if self.rest['Uhel'] == 30 else 1
        n_tr = math.ceil((math.pi * (DS-d_out)**2)/(4*heat_exchanger['t_t']**2*b1)) # ULOZIT PARAMETER
        A = (d_out * math.pi * n_tr) / (D1 ** 2 / 4) # ULOZIT PARAMETER 
        w_tube = self.tube.prutok / ((math.pi * heat_exchanger['d_in']**2) / 4 * self.tube.rho * n_tr) # ULOZIT PARAMETER 
        
        assert (self.tube.phase == 'liquid' or self.tube.phase == 'gas'), ''
        
        if self.tube.phase == 'liquid':
            if w_tube < 0.2 or w_tube > 1:
                raise Exception('')
        elif self.tube.phase == 'gas':
            if w_tube < 10 or w_tube > 30:
                raise Exception('')
        
        nu_tube = self.tube.eta / self.tube.rho
        Re_tube = w_tube * heat_exchanger['d_in'] / nu_tube
        Pr_tube = self.tube.eta * self.tube.cp / self.tube.lamb
        if self.tube.t1 > self.tube.t2:
            Nuss_tube = 0.023 * Re_tube ** 0.8 * Pr_tube ** 0.3
        else:
            Nuss_tube = 0.023 * Re_tube ** 0.8 * Pr_tube ** 0.4
        
        alpha_tube = Nuss_tube * self.tube.lamb / heat_exchanger['d_in']
        # ----------------------------------------------------------------

        # -----------------------MEZITRUBKOVY PROSTOR-----------------------
        x6 = t_t1 / d_out # 56
        x7 = t_t2 / d_out # 56
        if x7 < 1:
            gamma = 1 - math.pi /(4*x6*x7)
        else:
            gamma = 1 - math.pi /(4*x6) 
        
        fi_vs_a = 2 * math.acos((2/DS-d_out) * (h_p-D1/2)) # 71
        n_tr_v = round((DS - d_out)**2/(8*heat_exchanger['t_t']**2*b1)*(fi_vs_a - math.sin(fi_vs_a))) # 58 

        S_2N = (heat_exchanger['t_p'] - s_p)*D1 # 55 velikost nezaplneneho prostoru mezi prepazkami
        w_shell = self.shell.prutok / (S_2N * self.shell.rho * gamma) # ULOZIT PARAMETER
        
        assert (self.shell.phase != 'liquid' or self.shell.phase != 'gas'), ''
                
        if self.shell.phase == 'liquid':
            if w_shell < 0.2 or w_shell > 0.8:
                raise Exception('')
        elif self.shell.phase == 'gas':
            if w_shell < 5 or w_shell > 15:
                raise Exception('')
        
        l_char = math.pi * d_out / 2
        Re_shell = (w_shell * l_char)/(self.shell.eta / self.shell.rho)
        Pr_shell = self.shell.cp * self.shell.eta / self.shell.lamb
        Nuss_shell_lam = 0.664 * Re_shell ** 0.5 * Pr_shell ** (1/3)
        Nuss_shell_turb = (0.037*Re_shell**0.8*Pr_shell)/(1+2.443*Re_shell**(-0.1)*(Pr_shell**(2/3)-1)) # 55

        # KOREKCNI FAKTORY
        y2 = 1 # TODO KONTROLA pomer prestupu tepla medii -> podle toho stena trubky
        
        if self.rest['Uhel'] == 90: #55
            y3 = 1 + 0.7/gamma**1.5*(x7/x6-0.3)/(x7/x6+0.7)**2
        else:
            y3 = 1+2/(3*x7) # 56 prevod soucinitele prestupu tepla z rady na svazek trubek 2 druhy pro 90 stupnu a pro zbytek
        
        y4 = 1 # TODO KONTROLA

        x8 = n_tr_v / n_tr # 56

        assert(x8 < 0.8), ''
        
        y5 = 1 - x8 + 0.524 * x8**0.32 # 57 0.2 < t_p/D1 < 1; x8 < 0.8

        y6 = 1 # 58 TODO bud 1 mm nebo zjistit a zohlednit

        if self.rest['Uhel'] == 30 or self.rest['Uhel'] == 90:
            S_2Z = ((D1-DS)+(heat_exchanger['t_t'] - d_out)*((DS-d_out)/t_t1))*(heat_exchanger['t_p'] - s_p) # 58-59 
        elif self.rest['Uhel'] == 45 or self.rest['Uhel'] == 60:
            S_2Z = ((D1-DS)+(heat_exchanger['t_t'] - d_out)*((DS-d_out)/(t_t1/2)))*(heat_exchanger['t_p'] - s_p) # 58-59
        S_sS = (D1 - DS - (heat_exchanger['t_t'] - d_out))*(heat_exchanger['t_p'] - s_p)
        if Re_shell < 100:
            y7 = math.exp(0-1.5*S_sS/S_2Z) # 59 TODO zanedbani tesnicich list?
        else:
            y7 = math.exp(0-1.35*S_sS/S_2Z) # 59 TODO zanedbani tesnicich list?
        y8 = 1 # TODO KONTROLA

        Nuss_shell = (0.3 +(Nuss_shell_lam**2 + Nuss_shell_turb**2)**0.5)*y2*y3*y4*y5*y6*y7*y8 # realny plati pri stehlik 60
        alpha_shell = Nuss_shell * self.shell.lamb / l_char
        # ----------------------------------------------------------------

        k = math.pi/(1/(alpha_tube * heat_exchanger['d_in']) + 1/(2*self.rest['Lambda'])*math.log(heat_exchanger['d_in']/d_out) + 1/(alpha_shell * d_out))
        if self.tube.t1 > self.tube.t2:
            delta_t_ln = ((self.tube.t1 - self.shell.t2) - (self.tube.t2 - self.shell.t1)) / math.log((self.tube.t1 - self.shell.t2) / (self.tube.t2 - self.shell.t1))
        else:
            delta_t_ln = ((self.shell.t1 - self.tube.t2) - (self.shell.t2 - self.tube.t1)) / math.log((self.shell.t1 - self.tube.t2) / (self.shell.t2 - self.tube.t1))
        L_max = round(self.rest['Q'] / (k * delta_t_ln * n_tr), 3) # ULOZIT PARAMETER
        if L_max > 15 * D1 or L_max < 3 * D1 or L_max > self.rest['MaxDelka']: raise Exception('') # ERROR
        n_p = math.floor(L_max / heat_exchanger['t_p'] - 1) # Ulozit parameter

        # KONTROLA

        # TLAKOVE ZTRATY V TRUBKOVEM PROSTORU
        if Re_tube <= 2320:
            lamb_11 = 64/Re_tube
        else:
            x9 = (2.457 * math.log(1/(7/Re_tube)** 0.9 + 0.27 * 0.2 / (1000*heat_exchanger['d_in'])))**16 # absolutni drsnost ocele 0.2
            x10 = (37530/Re_tube)**16
            lamb_11 = 8*((8/Re_tube)**12 + 1/(x9 + x10)**1.5)**(1/12)
        
        z1 = L_max / heat_exchanger['d_in']
        z2 = 1 # TODO KONTROLA
        tube_delta_p = self.tube.rho * w_tube ** 2 / 2 * (0.7 + lamb_11 * z1 * z2)

        # TLAKOVE ZTRATY V MEZITRUBKOVEM PROSTORU
        lambda22 = self.pressure_drop(Re_shell, heat_exchanger['t_t'], d_out)
        w_2 = self.shell.prutok / (S_2Z * self.shell.rho)
        if Re_shell < 100:
            z3 = math.exp(0-4.5*S_sS/S_2Z) 
        else:
            z3 = math.exp(0-3.7*S_sS/S_2Z)
        z4 = 1 # TODO KONTROLA
        n_rp = (2 * h_p - D1)/t_t2
        shell_p_a = 2 * lambda22 * n_rp * (n_p - 1) * self.shell.rho * w_2 ** 2 * z2 * z3 * z4

        l_tn = 2*heat_exchanger['t_p'] - s_p
        if Re_shell < 100:
            z5 = 2*(2*heat_exchanger['t_p']/l_tn)**(2-1)
        else:
            z5 = 2*(2*heat_exchanger['t_p']/l_tn)**(2-0.2)
        n_rv = 0.8/t_t2 * ((D1 + (DS - d_out))/2 - h_p)
        shell_p_b = 2 * lambda22 * (n_rp + n_rv) * self.shell.rho * w_2 ** 2 * z2 * z3 * z5
        
        phi_vp = 2*math.acos(2*h_p/D1 -1)
        S_vN = math.pi* D1**2 / 4 *(phi_vp/(2*math.pi) - math.sin(phi_vp)/(2*math.pi))
        S_vZ = S_vN - n_tr_v*math.pi*d_out**2/4
        w_2v = self.shell.prutok / ((S_2Z * S_vZ)**0.5 * self.shell.rho)
        if Re_shell <= 100:
            d_hv = 4 * S_vZ/(n_tr_v*math.pi*d_out + math.pi * D1 * phi_vp/(2*math.pi))
            shell_p_c = n_p*(2*self.shell.rho*w_2v**2/2 +26 * self.shell.prutok *self.shell.eta/((S_vZ*S_vN)**0.5)*(n_rv/(heat_exchanger['t_t']-d_out)+heat_exchanger['t_p']/d_hv**2))*z4
        else:
            shell_p_c = n_p *((2 + 0.6 * n_rv)*self.shell.rho*w_2v**2/2)*z4
        shell_delta_p = shell_p_a + shell_p_b + shell_p_c
        delta_p = tube_delta_p + shell_delta_p # ULOZIT PARAMETER

        if delta_p > self.rest['MaxZtraty']:
            raise Exception('')
        
        # TODO pricit k objemu prepazky
        objem_vymeniku = math.pi*(n_tr*(d_out**2 - heat_exchanger['d_in']**2)+(heat_exchanger['shell']['D2']**2-D1**2))/4*L_max 

        result = {
            'vyska_prep': round(h_p,3),
            'tl_prep': s_p,
            'delka':L_max,
            'pocet_prepazek':n_p,
            'pocet_trubek': n_tr,
            'w_tube': round(w_tube,2),
            'w_shell': round(w_shell,2),
            'kompaktnost': round(A),
            'tlak_ztraty': round(delta_p),
            'hmotnost': objem_vymeniku * self.sizes.rho
        }
        return (dict(heat_exchanger), result)

    def baffle_thickness(self, shell_diameter, baffle_spacing):
        if shell_diameter < 0.377:
            if baffle_spacing < 0.3:
                return 0.004
            elif baffle_spacing < 0.45:
                return 0.005
            elif baffle_spacing < 0.6:
                return 0.006
            elif baffle_spacing < 0.8:
                return 0.01
            else: return 0.01
        elif shell_diameter < 0.7:
            if baffle_spacing < 0.3:
                return 0.005
            elif baffle_spacing < 0.45:
                return 0.006
            elif baffle_spacing < 0.6:
                return 0.01
            elif baffle_spacing < 0.8:
                return 0.01
            else: return 0.012
        elif shell_diameter < 1:
            if baffle_spacing < 0.3:
                return 0.006
            elif baffle_spacing < 0.45:
                return 0.008
            elif baffle_spacing < 0.6:
                return 0.01
            elif baffle_spacing < 0.8:
                return 0.012
            else: return 0.016
        else:
            if baffle_spacing < 0.3:
                return 0.006
            elif baffle_spacing < 0.45:
                return 0.01
            elif baffle_spacing < 0.6:
                return 0.012
            elif baffle_spacing < 0.8:
                return 0.016
            else: return 0.016
         
    def pressure_drop(self, reynolds, t_t, d_out):
        if self.rest['Uhel'] == 30 or self.rest['Uhel'] == 60:
            c2 = 7
            a2 = 0.5
            if reynolds > 10000:
                c1 = 0.372
                a1 = -0.123
            elif reynolds > 1000:
                c1 = 0.486
                a1 = -0.152
            elif reynolds > 100:
                c1 = 0.57
                a1 = -0.476
            elif reynolds > 10:
                c1 = 45.1
                a1 = -0.973
            elif reynolds > 0:
                c1 = 48
                a1 = -1
            else:
                raise Exception('')        
        elif self.rest['Uhel'] == 45:
            c2 = 6.59
            a2 = 0.52
            if reynolds > 10000:
                c1 = 0.303
                a1 = -0.126
            elif reynolds > 1000:
                c1 = 0.333
                a1 = -0.136
            elif reynolds > 100:
                c1 = 3.5
                a1 = -0.476
            elif reynolds > 10:
                c1 = 26.2
                a1 = -0.913
            elif reynolds > 0:
                c1 = 32
                a1 = -1
            else:
                raise Exception('')       
        elif self.rest['Uhel'] == 90:
            c2 = 6.3
            a2 = 0.378
            if reynolds > 10000:
                c1 = 0.391
                a1 = -0.148
            elif reynolds > 1000:
                c1 = 0.0815
                a1 = 0.022
            elif reynolds > 100:
                c1 = 6.09
                a1 = -0.602
            elif reynolds > 10:
                c1 = 32.1
                a1 = -0.963
            elif reynolds > 0:
                c1 = 35
                a1 = -1
            else:
                raise Exception('')       

        a = c2 / (1 + 0.14 * reynolds ** a2)
        lamb22 = c1 * (1.33/(t_t/d_out))**a * reynolds ** a1
        return lamb22

if __name__ == '__main__':
    trubka = {'M': 15/360, 'T1': 890 + 273.15, 'T2': 60.81 + 273.15, 'P': 100000, 'Rf': 10.0, 'Medium': [['H2O', 0.1],
    ['O2', 0.07],
    ['CO2', 0.12],
    ['Ar', 0.006],
    ['N2', 0.7]]}
    plast = {'M': 1.926, 'T1': 25 + 273.15, 'T2': 30 + 273.15, 'P': 400000, 'Rf': 10.0, 'Medium': [['H2O', 1]]}
    ostatni = {
        'Uhel': 30.0,
        'Lambda': 50.0,
        'MaxDelka': 6.0,
        'MaxSirka': 1.2,
        'MaxZtraty': 50000.0,
    }
    everything = Calculate(trubka, plast, ostatni)
    everything.calculate_all()