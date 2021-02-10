import math
import time
from CoolProp.CoolProp import PropsSI
from elements import Element
cas_startu = time.time()

uhel = 30
D_max = 0.7 # [m]
L_max = 5 # [m]
m_sp = 150/3600 # [kg/s]
t_sp_in = 890 +273.15 # [Kelvin]
p_sp = 10**5 # [Pascal]
lamb_ocel = 50 #ZVOLENO
# procentualni zastoupeni prvku ve spalinach
# poradi [O2, CO2, N2, Ar, H2O]s
x = [
    0.07,
    0.12,
    0.7,
    0.006,
    0.1
]

# molarni hmotnost prvku
# poradi [O2, CO2, N2, Ar, H2O]
M = [
    32,
    44,
    28,
    40,
    18
]

Spaliny = Element(x, M)
w = Spaliny.w()
# voda
t_w_in = 25 + 273.15 # [Celsius]
t_w_out = 30 +273.15 # [Celsius]
t_w_str = (t_w_in+t_w_out)/2
p_w = 4*10**5 # [Pascal]

# teplota rosneho bodu
pp_sp_H2O = p_sp * x[4]
t_sat = PropsSI('T','P',pp_sp_H2O,'Q',0,'Water')
t_sp_out = t_sat + 15

# entalpie spalin na vstupu[J/kg]
i_sp_in = Spaliny.enthalpy(t_sp_in, p_sp)

# entalpie spalin na vystupu [J/kg]
i_sp_out = Spaliny.enthalpy(t_sp_out, p_sp)

# entalpie vody
i_w_in = Spaliny.enthalpy_w(t_w_in, p_w)
i_w_out = Spaliny.enthalpy_w(t_w_out, p_w)

# Hmotnostni prutok vody
Q_sp = m_sp*(i_sp_in - i_sp_out) #[W]
m_w = Q_sp / (i_w_out - i_w_in) # [kg/s]
Qv = 1.05*Q_sp #ZVOLENO

# stredni teplota spalin
t_sp_str = (t_sp_in + t_sp_out)/2

slozky = [
    ['O2', 0],
    ['CO2', 0],
    ['N2', 0],
    ['Ar', 0],
    ['H2O', 0]
]

# hustota spalin
rho = 0 # [kg/m**3]
for i in slozky:
    i[1] = PropsSI('D', 'T', t_sp_str, 'P', p_sp, i[0])
for i in range(5):
    rho += slozky[i][1] * x[i]

#dynamicka viskozita spalin
eta = 0 #[m**2/s]
for i in slozky:
    i[1] = PropsSI('V', 'T', t_sp_str, 'P', p_sp, i[0])
for i in range(5):
    eta += slozky[i][1] * x[i]

# merna tepelna kapacita spalin
cp = 0 #[W/kgK]
for i in slozky:
    i[1] = PropsSI('CPMASS', 'T', t_sp_str, 'P', p_sp, i[0])
for i in range(5):
    cp += slozky[i][1] * x[i]

#soucinitel tepelne vodivosti
lamb = 0 # [W/mK]
for i in slozky:
    i[1] = PropsSI('L', 'T', t_sp_str, 'P', p_sp, i[0])
for i in range(5):
    lamb += (slozky[i][1] * x[i] / M[i]**0.5)/(x[i]/M[i]**0.5)
lamb *= 0.2
lamb2 = 0

for i in range(5):
    lamb2 += (slozky[i][1])
lamb2 *= 0.2
rho_w = PropsSI('D', 'T', t_w_str, 'P', p_w, 'Water')
eta_w = PropsSI('V', 'T', t_w_str, 'P', p_w, 'Water')
lamb_w = PropsSI('L', 'T', t_w_str, 'P', p_w, 'Water')
cp_w = PropsSI('CPMASS', 'T', t_w_str, 'P', p_w, 'Water')
# ---------------------------ROZMERY---------------------------
# Geometricky navrh ZVOLENO 
# trubky Stehlik 123
'''
d_in = 0.01 #[m] minimalne 12 mm kvuli zanaseni, maximalne 30 kvuli kompaktnosti PARAMETR
tl = 0.001 #[m] 1,5-5 mm  PARAMETR cim vic tim vetsi tlakova ztrata, ostatni zmeny male
uhel = 30 # [deg]
# plast Stehlik 123
D1 = 0.2 #[m] vnitrni prumer TAB
DS = 0.94 * D1 # [m] prumer trubkoveho svazku TAB
t_p = 1.5 * D1 #[m] vzdalenost mezi prepazkami ZVOLENO PARAMETR
d_out = d_in + 2*tl #[m]
t_t = 1.5*d_out #[m] roztec trubek 1,25-1,33-1,5 * d_in ZVOLENO PARAMETR alespon 6mm
'''
tlakova_ztrata = 100000
Pocet = 0
for D1 in [1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2]:
    pocet = 0
    if D1 > D_max: continue
    DS = D1 - 0.016
    for d_in in [0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028, 0.03]:
        for tl in [0.001, 0.002, 0.003, 0.004, 0.005]:
            d_out = d_in + 2 * tl
            for t_p in [0.2, 0.4, 0.6, 0.8, 1]:
                t_p *= D1
                for t_t in [1.25, 1.33, 1.5]:
                    # stehlik 125
                    h_p = 0.75*DS # 50-75% * DS zakryty prostor prepazkou PARAMETR
                    s_p = 0.004 #[m] tloustka prepazek TAB
                    t_t *= d_out
                    # stehlik 56
                    if uhel == 30:
                        c1 = 1
                        c2 = 0.866
                    elif uhel == 45:
                        c1 = 1.414
                        c2 = 0.707
                    elif uhel == 60:
                        c1 = 1.732
                        c2 = 0.5
                    elif uhel == 90:
                        c1 = 1
                        c2 = 1
                    else:
                        print('Error spatny uhel')
                    t_t1 = c1*t_t #[m] pricna roztec trubek TAB ruzne druhy pro ruzne stupne viz stehlik 56
                    t_t2 = c2*t_t #[m] podelna roztec trubek TAB ruzne druhy pro ruzne stupne viz stehlik 56
                    # pocet trubek Stehlik 57,58,71
                    if uhel == 30:
                        b1 = 0.866
                    else:
                        b1 = 1 
                    n_tr = math.ceil((math.pi * (DS-d_out)**2)/(4*t_t**2*b1)) # 0.866 = 30 stupnu, 1 = 45,60,90 stupnu
                    
                    
                    # -----------------------TRUBKOVY PROSTOR-----------------------
                    w_sp = m_sp / ((math.pi*d_in**2)/4*rho*n_tr) #rychlost spalin CHYBI POCET CHODU 5-10
                    if w_sp > 15 or w_sp < 8: continue
                    nu_sp = eta / rho # kinematicka viskozita
                    Re_sp = (w_sp*d_in)/nu_sp # Reynolds
                    Pr_sp = eta * cp / lamb # Prandlt
                    if t_sp_in > t_sp_out:
                        n = 0.3
                    else:
                        n = 0.4    
                    Nuss_sp = 0.023 * Re_sp**0.8 * Pr_sp**n # Nusselt

                    alpha_sp = Nuss_sp * lamb / d_in


                    fi_vs_a = 2 * math.acos((2/DS-d_out) * (h_p-D1/2)) # 71
                    n_tr_v = round((DS - d_out)**2/(8*t_t**2*b1)*(fi_vs_a - math.sin(fi_vs_a))) # 58 

                    S_2N = (t_p - s_p)*D1 # 55 velikost nezaplneneho prostoru mezi prepazkami
                    l_char = math.pi * d_out / 2 # 55 charakteristicky rozmer
                    s_tt = t_t - d_out
                    S_sS = (D1 - DS - s_tt)*(t_p - s_p)
                    if uhel == 30 or uhel == 90:
                        S_2Z = ((D1-DS)+s_tt*((DS-d_out)/t_t1))*(t_p-s_p) # 58-59 zalezi na stupnich 
                    elif uhel == 45 or uhel == 60:
                        S_2Z = ((D1-DS)+s_tt*((DS-d_out)/(t_t1/2)))*(t_p-s_p) # 58-59 zalezi na stupnich
                    n_rp = (2*h_p - D1)/t_t2 # 57
                    n_rv = (0.8/t_t2)*((D1 + DS - d_out)/2-h_p) # 70
                    fi_vp = 2*math.acos(2*h_p/D1-1) # 58
                    l_tn = 2*t_p-s_p

                    # podelne proudeni 71-72
                    S_vn = (math.pi*D1**2)/4*(fi_vp/(2*math.pi)-math.sin(fi_vp)/(2*math.pi)) # UPRAVENO KVULI RADIANUM
                    S_vz = S_vn-n_tr_v*((math.pi*d_out**2)/4)

                    # -----------------------MEZITRUBKOVY PROSTOR-----------------------
                    # Stehlik
                    x6 = t_t1 / d_out # 56
                    x7 = t_t2 / d_out # 56
                    if x7 < 1:
                        gamma = 1 - math.pi /(4*x6*x7)
                    else:
                        gamma = 1 - math.pi /(4*x6) 
                    x8 = n_tr_v / n_tr # 56

                    w_w = m_w / (S_2N * rho_w * gamma) # 55 rychlost vody, vztah zalezi na x7
                    if w_w > 0.8 or w_w < 0.2: continue
                    Re_w = (w_w*l_char)/(eta_w/rho_w) # 55
                    Pr_w = cp_w*eta_w/lamb_w # 55
                    Nu_lam = 0.664*Re_w**0.5*Pr_w**(1/3) # 55
                    Nu_turb = (0.037*Re_w**0.8*Pr_w)/(1+2.443*Re_w**(-0.1)*(Pr_w**(2/3)-1)) # 55
                    # stehlik 56-60
                    y2 = 1 # zmena latkovych vlastnosti v mezni vrstve DOPLNIT TODO

                    if uhel == 90: #55
                        y3 = 1 + 0.7/gamma**1.5*(x7/x6-0.3)/(x7/x6+0.7)**2
                    else:
                        y3 = 1+2/(3*x7) # 56 prevod soucinitele prestupu tepla z rady na svazek trubek 2 druhy pro 90 stupnu a pro zbytek
                    if Re_w < 100: # 57 zalezi na Re (3 vztahy)
                        n_rc = n_rp *(3-1) # TODO 3 by mela by n_p (pocet prepazek)
                        y4 = 1.51/n_rc**0.18
                        if Re_w > 20:
                            y4 = y4 + (20-Re_w)/80*(y4 - 1)
                    else:
                        y4 = 1

                    y5 = 1 - x8 + 0.524 * x8**0.32 # 57 0.2 < t_p/D1 < 1; x8 < 0.8
                    if x8 > 0.8 or t_p/D1 < 0.2 or t_p/D1 > 1:
                        print('ERROR x8')

                    y6 = 1 # 58 TODO zjistit dulezitost, spatne by se zjistovalo
                    if Re_w < 100:
                        y71 = 1.5
                    else:
                        y71 = 1.35
                    y7 = math.exp(0-y71*S_sS/S_2Z) # 59 TODO zanedbani tesnicich list?
                    if S_sS/S_2Z > 0.5:
                        print('ERROR y7')
                    if Re_w <= 100:
                        y81 = 0.33
                    else:
                        y81 = 0.6
                    y8 = 1 # 59 TODO zjistit dulezitost pocet prepazek

                    Nuss_w = (0.3 +(Nu_lam**2 + Nu_turb**2)**0.5)*y2*y3*y4*y5*y6*y7*y8 # realny plati pri stehlik 60
                    # TODO
                    # plati pokud n_rc = (n_p-1)*n_rp > 10
                    # Re (10;10**6)
                    # Pr (0.6; 10**3)
                    alpha_w = Nuss_w * lamb_w / l_char
                    # TODO zapocitani zanaseni
                    k = math.pi/(1/(alpha_sp * d_in) + 1/(2*lamb_ocel)*math.log(d_in/d_out) + 1/(alpha_w * d_out)) # prostup tepla
                    # TODO korekcni faktor pro vice chodu
                    delta_t_ln = ((t_sp_in-t_w_out)-(t_sp_out-t_w_in))/math.log((t_sp_in-t_w_out)/(t_sp_out-t_w_in)) # teplotni dif
                    l_max = round(Qv/(k*n_tr*delta_t_ln),3) # delka vymeniku minimalne 3*D1
                    if l_max > 15*D1 or l_max < 3*D1: continue
                    A = l_max * d_in**2 * math.pi /4 * n_tr
                    Qv_1 = l_max*k*n_tr*delta_t_ln
                    n_p = math.floor(l_max/t_p -1) # pocet prepazek

                    parametry = [D1, d_in, tl, t_p, t_t]
                    # TLAKOVE ZTRATY trubkovy prostor 66-67
                    # nezapocitani vstupniho a vystupniho hrdla
                    # Prepazkovy prostor trenim
                    lamb11 = 64 / Re_sp #66 zalezi na Re
                    # TODO turbulentni proudeni
                    '''
                    if Re_sp > 2320:
                        print('ERROR Turbulentni proudeni lamb11')
                    '''
                    z_1 =l_max / d_in
                    z_2 = 1 # 67 Re mezni vrstva TODO zjistit dulezitost
                    delta_p_t = lamb11 * rho * w_sp**2 / 2 * z_1 * z_2
                    # Prepazkovy prostor mistni
                    delta_p_m = 0.7 * rho * w_sp ** 2 / 2 # 67 TODO POCET CHODU
                    # Celkova ztrata v trubkovem prostoru
                    delta_p = delta_p_m + delta_p_t

                    # TLAKOVE ZTRATY mezitrubkovy prostor 68-71
                    # VELIKEJ PROBLEM URCIT KONSTANTY, TABULKA JAK KRAVA NA STRANE 69 TODO
                    c1 = 0.486
                    c2 = 7
                    a1=-0.152
                    a2=0.5
                    w_zw = m_w / (S_2Z*rho_w)
                    Re_zw = (w_zw*l_char)/(eta_w/rho_w)
                    a=c2/(1+0.14*Re_zw**a2)

                    lamb22 = c1*(1.33/(t_t/d_out))**a*Re_zw**a1

                    z_2w = 1 # 69 TODO
                    if Re_w < 100:
                        z_31w = 4.5
                    else:
                        z_31w = 3.7
                    z_3w = math.exp(0-z_31w*(S_sS/S_2Z)) # 70
                    z_4w = 1 # 70 TODO

                    delta_p_to = 2 * lamb22 * n_rp * (n_p - 1)*rho_w *w_zw**2 * z_2w* z_3w*z_4w
                    # Neprepazkovany prostor
                    if Re_w < 100:
                        z_51w = 1
                    else:
                        z_51w = 0.2
                    z_5w = 2*(2*t_p/l_tn)**(2-z_51w) # 70
                    delta_p_tn = 2*lamb22*(n_rp + n_rv)* rho_w*w_zw**2* z_2w* z_3w* z_5w
                    w_wv = m_w/((S_2Z* S_vz)**0.5*rho_w)
                    delta_p_tv = n_p*((2+0.6*n_rv)*rho_w*w_wv**2)/2*z_4w # TODO laminarni proudeni
                    delta_p_w = delta_p_to + delta_p_tn + delta_p_tv

                    pocet += 1
                    Pocet += 1
                    if delta_p_w + delta_p < tlakova_ztrata:
                        min_ztrata = parametry
                        tlakova_ztrata = delta_p + delta_p_w
    print('D1: {} a pocet: {}'.format(D1, pocet))
print('Celkove: {}'.format(Pocet))
print(min_ztrata)


print(time.time()-cas_startu)