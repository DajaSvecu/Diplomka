import CoolProp
import math
from elements import Element
from CoolProp.CoolProp import PropsSI
m_sp = 150/3600 # [kg/s]
t_sp_in = 890 +273.15 # [Kelvin]
p_sp = 10**5 # [Pascal]
# procentualni zastoupeni prvku ve spalinach
# poradi [O2, CO2, N2, Ar, H2O]
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

rho_w = PropsSI('D', 'T', t_w_str, 'P', p_w, 'Water')
eta_w = PropsSI('V', 'T', t_w_str, 'P', p_w, 'Water')
lamb_w = PropsSI('L', 'T', t_w_str, 'P', p_w, 'Water')
cp_w = PropsSI('CPMASS', 'T', t_w_str, 'P', p_w, 'Water')

# -----------------------TRUBKOVY PROSTOR-----------------------
# Geometricky navrh ZVOLENO 
# trubky Stehlik 123
d_in = 0.015 #[m] minimalne 12 mm kvuli zanaseni, maximalne 30 kvuli kompaktnosti PARAMETR
tl = 0.002 #[m] 1,5-5 mm  PARAMETR
d_out = d_in + 2*tl #[m]

# plast Stehlik 123
D1 = 0.2 #[m] vnitrni prumer TAB
DS = 0.186 # [m] prumer trubkoveho svazku TAB
t_p = 1*D1 #[m] vzdalenost mezi prepazkami ZVOLENO PARAMETR
t_t = 1.5*d_in #[m] roztec trubek 1,25-1,5 * d_in ZVOLENO PARAMETR
# stehlik 56
t_t1 = t_t #[m] pricna roztec trubek TAB
t_t2 = 0.866*t_t #[m] podelna roztec trubek TAB
# stehlik 125
s_p = 0.004 #[m] tloustka prepazek TAB

# pocet trubek Stehlik 57,58,71
n_tr = round((math.pi * (DS-d_out)**2)/(4*t_t**2*0.866)) # 0.866 = 30 stupnu, 1 = 45,60,90 stupnu
h_p = 0.5*DS # 50-75% * DS zakryty prostor prepazkou PARAMETR
fi_vs_a = 2 * math.acos((2/DS-d_out) * (h_p-D1/2)) # 71
n_tr_v = round((DS - d_out)**2/(8*t_t**2*0.866)*(fi_vs_a - math.sin(fi_vs_a))) # 58 

# rychlost v trubkach
w_sp = m_sp / ((math.pi*d_in**2)/4*rho*n_tr) #rychlost spalin CHYBI POCET CHODU 5-10
nu_sp = eta / rho # kinematicka viskozita
Re_sp = (w_sp*d_in)/nu_sp # Reynolds
Pr_sp = eta * cp / lamb # Prandlt

Nuss_sp = 0.023 * Re_sp**0.8 * Pr_sp**0.4 # Nusselt

alpha_sp = Nuss_sp * lamb / d_in

# -----------------------MEZITRUBKOVY PROSTOR-----------------------
# Stehlik
S_n = (t_p - s_p)*D1 # 55 velikost nezaplneneho prostoru mezi prepazkami
l_char = math.pi * d_out / 2 # 55 charakteristicky rozmer
x6 = t_t1 / d_out # 56
x7 = t_t2 / d_out # 56
x8 = n_tr_v / n_tr # 56
w_w = m_w / (S_n * rho_w * (1-math.pi/(4*x6))) # 55 rychlost vody 0,2-0,8

Re_w = (w_w*l_char)/(eta_w/rho_w) # 55
Pr_w = cp_w*eta_w/lamb_w # 55

# stehlik 56-60
y2 = 1 # zmena latkovych vlastnosti v mezni vrstve
y3 = 1+2/(3*x7) # prevod soucinitele prestupu tepla z rady na svazek trubek
y4 = 1 # nepriznivy tvar teploniho profilu v proudu pracovni latky pri laminarnim proudeni
y5 = 1 - x8 + 0.524 * x8**0.32 # 0.2 < t_p/D1 < 1; x8 < 0.8 NEPLATI ZEPTAT SE
y6 = 1 # zanedbavame
s_tt = t_t - d_out
S_sS = (D1 - DS - s_tt)*(t_p - s_p)
S_2Z = (D1-DS)+(t_t-d_out)*((DS-d_out)/t_t1)*(t_p-s_p)
y7 = math.exp(0-1.35*S_sS/S_2Z) # nemame tesnici listy
y8 = 1 # vliv neoprepazkovany prostor

Nu_lam = 0.664*Re_w**0.5*Pr_w**(1/3) # 55
Nu_turb = (0.037*Re_w**0.8*Pr_w)/(1+2.443*Re_w**(-0.1)*(Pr_w**(2/3)-1)) # 55

# stehlik 60
Nuss_w = (0.3 +(Nu_lam**2 + Nu_turb**2)**0.5)*y2*y3*y4*y5*y6*y7*y8 # realny
Nuss_w_i = (0.3 +(Nu_lam**2 + Nu_turb**2)**0.5)*y2*y3*y4 # idealni
ucinnost = Nuss_w/Nuss_w_i # blizeni se k idealnimu
alpha_w = Nuss_w * lamb_w / l_char
lamb_ocel = 50 #ZVOLENO
k = math.pi/(1/(alpha_sp * d_in) + 1/(2*lamb_ocel)*math.log(d_in/d_out) + 1/(alpha_w * d_out)) # prostup tepla
delta_t_ln = ((t_sp_in-t_w_out)-(t_sp_out-t_w_in))/math.log((t_sp_in-t_w_out)/(t_sp_out-t_w_in)) # teplotni dif
Qv = Q_sp*1.05 #ZVOLENO
l_max = Qv/(k*n_tr*delta_t_ln) # delka vymeniku minimalne 3*D1
n_p = l_max/t_p -1 # pocet prepazek

# TLAKOVE ZTRATY trubkovy prostor 66-67
# Prepazkovy prostor trenim
lamb11 = 64 / Re_sp
z_1 =l_max / d_in
z_2 = 1 # zanedbavame
delta_p_t = lamb11 * rho * w_sp**2 / 2 * z_1 * z_2
# Prepazkovy prostor mistni
delta_p_m = 0.7 * rho * w_sp ** 2 / 2
# Celkova ztrata v trubkovem prostoru
delta_p = delta_p_m + delta_p_t
# TLAKOVE ZTRATY mezitrubkovy prostor 68-71
c1 = 0.486
c2 = 7
a1=-0.152
a2=0.5
w_zw = m_w / (S_2Z*rho_w)
Re_zw = (w_zw*l_char)/(eta_w/rho_w)
a=c2/(1+0.14*Re_zw**a2)

lamb22 = c1*(1.33/(t_t/d_out))**a*Re_zw**a1
 
z_2w = 1
z_3w = math.exp(-3.7*(S_sS/S_2Z))
z_4w = 1
n_rp = (2*h_p - D1)/t_t2

delta_p_to = 2 * lamb22 * n_rp * (n_p - 1)*rho_w *w_zw**2 * z_2w* z_3w*z_4w
# Neprepazkovany prostor
n_rv = (0.8/t_t2)*((D1 +DS - d_out)/2-h_p)
z_5w = 2*(2*t_p/(2*t_p-s_p))**(2-0.2)
delta_p_tn = 2*lamb22*(n_rp + n_rv)* rho_w*w_zw**2* z_2w* z_3w* z_5w
# podelne proudeni
fi_vp = 2*math.acos(2*h_p/D1-1)

s_vn = (math.pi*D1**2)/4*(fi_vp/(2*math.pi)-math.sin(fi_vp)/(2*math.pi)) # UPRAVENO KVULI RADIANUM
s_vz = s_vn-n_tr_v*((math.pi*d_out**2)/4)
w_wv = m_w/((S_2Z* s_vz)**0.5*rho_w)
delta_p_tv = n_p*((2+0.6*n_rv)*rho_w*w_wv**2)/2*z_4w
 
delta_p_w = delta_p_to + delta_p_tn + delta_p_tv

Dulezite = {
    'roztec': t_t,
    'Qv': Qv,
    'deltaPw': delta_p_w,
    'deltaPsp': delta_p,
    'D1':D1,
    'lmax': l_max,    
    'n_tr': n_tr,
    'w_sp': w_sp,
    'w_w': w_w,
    'pocet prepazek': n_p,
    'alpha_w':alpha_w,
    'alpha_sp':alpha_sp    
}