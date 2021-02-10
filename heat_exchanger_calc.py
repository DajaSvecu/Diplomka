from Medium import Medium
from heat_exchanger_inputs import Sizes
import math

class Calculate():
    def __init__(self, tube:dict, shell:dict, rest:dict):
        self.tube = Medium(tube)
        self.shell = Medium(shell)
        self.rest = rest
        self.sizes = Sizes()
        if self.tube.q > 0.99 * self.shell.q and self.tube.q < 1.01 * self.shell.q:
            self.rest['Q'] = 1.08 * (self.tube.q + self.shell.q) / 2
        #else: # ERROR
    
    def calculate_all(self) -> None:
        heat_exchanger = {}
        for shell in self.sizes.list_of_shells:
            if shell['D2'] > self.rest['MaxSirka']: continue # podminka pro maximalni sirku
            heat_exchanger['shell'] = shell
            for tube in self.sizes.list_of_tubes:
                heat_exchanger['d_in'] = tube['d_in']
                for tl in tube['wall']:
                    heat_exchanger['tl'] = tl
                    for t_p in {0.2, 0.4, 0.6, 0.8, 1}:
                        heat_exchanger['t_p'] = t_p
                        for t_t in {1.25, 1.33, 1.5}:
                            heat_exchanger['t_t'] = t_t
                            try:
                                heat_exchanger_result = self.calculate_heat_exchanger(heat_exchanger)
                            except AssertionError:
                                pass
                            except Exception:
                                pass
                            else:
                                print('This one work out!')

                            

    def calculate_heat_exchanger(self, heat_exchanger:dict) -> dict:
        D1 = heat_exchanger['shell']['D2'] - 2 * heat_exchanger['shell']['wall']
        DS = D1 - 0.016
        d_out = heat_exchanger['d_in'] + 2 * heat_exchanger['tl']
        heat_exchanger['t_p'] *= D1
        heat_exchanger['t_t'] *= d_out
        h_p = 0.75 * DS # ULOZIT PARAMETER
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
            if w_tube < 0.3 or w_tube > 1:
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
        y2 = 1 # zanedbavame
        
        if self.rest['Uhel'] == 90: #55
            y3 = 1 + 0.7/gamma**1.5*(x7/x6-0.3)/(x7/x6+0.7)**2
        else:
            y3 = 1+2/(3*x7) # 56 prevod soucinitele prestupu tepla z rady na svazek trubek 2 druhy pro 90 stupnu a pro zbytek
        
        y4 = 1 # KONTROLA

        x8 = n_tr_v / n_tr # 56

        assert(x8 < 0.8), ''
        
        y5 = 1 - x8 + 0.524 * x8**0.32 # 57 0.2 < t_p/D1 < 1; x8 < 0.8

        y6 = 1 # 58 TODO zjistit dulezitost, spatne by se zjistovalo

        if self.rest['Uhel'] == 30 or self.rest['Uhel'] == 90:
            S_2Z = ((D1-DS)+(heat_exchanger['t_t'] - d_out)*((DS-d_out)/t_t1))*(heat_exchanger['t_p'] - s_p) # 58-59 
        elif self.rest['Uhel'] == 45 or self.rest['Uhel'] == 60:
            S_2Z = ((D1-DS)+(heat_exchanger['t_t'] - d_out)*((DS-d_out)/(t_t1/2)))*(heat_exchanger['t_p'] - s_p) # 58-59
        S_sS = (D1 - DS - (heat_exchanger['t_t'] - d_out))*(heat_exchanger['t_p'] - s_p)
        if Re_shell < 100:
            y7 = math.exp(0-1.5*S_sS/S_2Z) # 59 TODO zanedbani tesnicich list?
        else:
            y7 = math.exp(0-1.35*S_sS/S_2Z) # 59 TODO zanedbani tesnicich list?
        y8 = 1 # KONTROLA

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

        #kontrola

        # TLAKOVE ZTRATY V TRUBKOVEM PROSTORU
        if Re_tube <= 2320:
            lamb_11 = 64/Re_tube
        else:
            x9 = (2.457 * math.log(1/(7/Re_tube)** 0.9 + 0.27 * 0.2 / (1000*heat_exchanger['d_in'])))**16 # absolutni drsnost ocele 0.2
            x10 = (37530/Re_tube)**16
            lamb_11 = 8*((8/Re_tube)**12 + 1/(x9 + x10)**1.5)**(1/12)
        
        z1 = L_max / heat_exchanger['d_in']
        z2 = 1 # zanedbavame
        tube_delta_p = self.tube.rho * w_tube ** 2 / 2 * (0.7 + lamb_11 * z1 * z2)

        # TLAKOVE ZTRATY V MEZITRUBKOVEM PROSTORU

        


        return 0

        
        
        
        
        
        
        



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
         




if __name__ == '__main__':
    trubka = {'M': 0.1, 'T1': 623.15, 'T2': 523.15, 'P': 100000.0, 'Rf': 10.0, 'Medium': [['Water', 1]]}
    plast = {'M': 0.48, 'T1': 298.15, 'T2': 308.15, 'P': 100000.0, 'Rf': 10.0, 'Medium': [['Water', 1]]}
    ostatni = {
        'Uhel': 30.0,
        'Lambda': 47.0,
        'MaxDelka': 5.0,
        'MaxSirka': 1.2,
        'MaxZtraty': 50000.0,
    }
    everything = Calculate(trubka, plast, ostatni)
    everything.calculate_all()