from heat_exchanger_calc import Calculate
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import math

class Calculate_param_study():
    def __init__(self, tube, shell, rest):
        self.basic = Calculate(tube, shell, rest)
        self.results = {}
        self.params = {
            't_p' : np.arange(0.2,0.95,0.05),
            'h_p' : np.arange(0.55,0.925,0.025),
            'tl' : np.arange(0.05,0.26,0.01),
            #'tl' : np.arange(0.5,2.6,0.1),
            't_t' : np.arange(1.2,2,0.05),
            'length' : np.arange(1,6,0.25),
            'd_out' : self.basic.sizes.list_of_tubes,
            'shell' : self.basic.sizes.list_of_shells,
        }
        self.param_names = {
            't_p' : 'Roztec prepazky [mm]',
            'h_p' : 'Vyska prepazky [mm]',
            'tl' : 'Tloustka steny trubky [mm]',
            't_t' : 'Roztec trubek [mm]',
            'length' : 'Delka trubek [m]',
            'd_out' : 'Vnejsi prumer trubek [mm]',
            'shell' : 'DN [-]', 
        }
        self.param_study_t_p()
        self.param_study_h_p()
        self.param_study_tl()
        self.param_study_t_t()
        self.param_study_shell()
        self.param_study_length()     
        self.param_study_d_out()



    def param_study_t_p(self):
        vysledky = []
        for t_p in self.params['t_p']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': 0.0127,
            'tl': 0.000889,
            'length': 1.2192,
            't_t': 0.015875,
            't_p':0.202642*t_p,
            'h_p':0.202642*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['t_p'] = vysledky
    def param_study_h_p(self):
        vysledky = []
        for h_p in self.params['h_p']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': 0.0127,
            'tl': 0.000889,
            'length': 1.2192,
            't_t': 0.015875,
            't_p':0.202642*0.5,
            'h_p':0.202642*h_p
            }, False)
            vysledky.append(vysledek)
        self.results['h_p'] = vysledky
    def param_study_tl(self):
        vysledky = []
        for tl in self.params['tl']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': 0.0127,
            'tl': tl*0.0127,
            'length': 1.2192,
            't_t': 0.019, # MEGA ZAJIMAVE, PRI MALE ROZTECI NEMA TLOUSTKA TEMER VLIV, PRI VETSI ZACINA KLESAT
            't_p':0.202642*0.5,
            'h_p':0.202642*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['tl'] = vysledky
    def param_study_t_t(self):
        vysledky = []
        for t_t in self.params['t_t']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': 0.0127,
            'tl': 0.000889,
            'length': 1.2192,
            't_t': 0.0127*t_t,
            't_p':0.202642*0.5,
            'h_p':0.202642*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['t_t'] = vysledky
    def param_study_length(self):
        vysledky = []
        for length in self.params['length']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': 0.0127,
            'tl': 0.000889,
            'length': length,
            't_t': 0.015875,
            't_p':0.202642*0.5,
            'h_p':0.202642*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['length'] = vysledky
    def param_study_d_out(self):
        vysledky = []
        for tube in self.params['d_out']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': self.basic.sizes.list_of_shells[0],
            'd_out': tube['d_out'],
            'tl': 0.000889,
            'length': 1.2192,
            't_t': tube['d_out']*1.25,
            't_p':0.202642*0.5,
            'h_p':0.202642*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['d_out'] = vysledky
    def param_study_shell(self):
        vysledky = []
        for shell in self.params['shell']:
            vysledek = self.basic.calculate_heat_exchanger({
            'shell': {
                'DN': shell['DN'], # [-]
                'wall': 0.008179, # [m]
                'D2': shell['D2'], # [m]
            },
            'd_out': 0.0127,
            'tl': 0.000889,
            'length': 1.2192,
            't_t': 0.015875,
            't_p':(shell['D2']-2*0.008179)*0.5,
            'h_p':(shell['D2']-2*0.008179)*0.8
            }, False)
            vysledky.append(vysledek)
        self.results['shell'] = vysledky


    def param_study_graph(self, param_name, study):
        fig, axs = plt.subplots(1,2, figsize=(15,10))
        data = {
            'x1' : [],
            'y1' : [],
            'x2' : [],
            'y2' : []
        }
        for result in self.results[param_name]:
            data['x1'].append(result[0][param_name])
            data['x2'].append(result[0][param_name])
            if study == 'vykon':
                data['y1'].append(result[1]['vykon'])    
            elif study == 'plocha':
                plocha = result[0]['d_out']*math.pi*result[1]['pocet_trubek']*result[0]['length']
                data['y1'].append(plocha)
            data['y2'].append(result[1]['tlak_ztraty'])
        
        if param_name != 'length':
            for i in range(len(data['x1'])):
                data['x1'][i] *=1000 
                data['x2'][i] *=1000
        axs[0].scatter(data['x1'], data['y1'])
        axs[1].scatter(data['x2'], data['y2'])

        z1 = np.polyfit(data['x1'], data['y1'],6)
        p1 = np.poly1d(z1)
        axs[0].plot(data['x1'],p1(data['x1']),'b--')
        axs[0].grid()
        z2 = np.polyfit(data['x2'], data['y2'],6)
        p2 = np.poly1d(z2)
        axs[1].plot(data['x2'],p2(data['x2']),'b--')
        axs[1].grid()

        if study == 'vykon':
            axs[0].set_title('Vykon')
        elif study =='plocha':
            axs[0].set_title('Teplosmenna plocha')
        axs[1].set_title('Tlakove ztraty')
        axs[0].set_xlabel(self.param_names[param_name]) 
        axs[1].set_xlabel(self.param_names[param_name]) 
        #axs[0].set_ylim(bottom=0)
        #axs[1].set_ylim(bottom=0)
        axs[0].yaxis.set_major_formatter(PercentFormatter(xmax=data['y1'][0], decimals=0))
        #axs[0].xaxis.set_major_formatter(PercentFormatter(xmax=202.642, decimals=0)) #202.642
        axs[1].yaxis.set_major_formatter(PercentFormatter(xmax=data['y2'][0], decimals=0))
        #axs[1].xaxis.set_major_formatter(PercentFormatter(xmax=202.642, decimals=0))
        plt.savefig("diplomka/pictures/{}_{}.png".format(param_name, study))
        plt.show()
    def param_study_shell_graph(self, study):
        fig, axs = plt.subplots(1, 2, figsize=(15,10))
        data = {
            'x1' : [],
            'y1' : [],
            'x2' : [],
            'y2' : []
        }
        for result in self.results['shell']:
            data['x1'].append(result[0]['shell']['DN'])
            data['x2'].append(result[0]['shell']['DN'])
            if study == 'vykon':
                data['y1'].append(result[1]['vykon'])    
            elif study == 'plocha':
                plocha = result[0]['d_out']*math.pi*result[1]['pocet_trubek']*result[0]['length']
                data['y1'].append(plocha)
            data['y2'].append(result[1]['tlak_ztraty'])

        axs[0].scatter(data['x1'], data['y1'])
        axs[1].scatter(data['x2'], data['y2'])
        z1 = np.polyfit(data['x1'], data['y1'],6)
        p1 = np.poly1d(z1)
        axs[0].plot(data['x1'],p1(data['x1']),'b--')
        z2 = np.polyfit(data['x2'], data['y2'],6)
        p2 = np.poly1d(z2)
        axs[1].plot(data['x2'],p2(data['x2']),'b--')

        if study == 'vykon':
            axs[0].set_title('Vykon')
        elif study =='plocha':
            axs[0].set_title('Teplosmenna plocha')
        axs[1].set_title('Tlakove ztraty')
        axs[0].set_xlabel(self.param_names['shell']) 
        axs[1].set_xlabel(self.param_names['shell']) 
        axs[0].yaxis.set_major_formatter(PercentFormatter(xmax=data['y1'][0], decimals=0))
        axs[1].yaxis.set_major_formatter(PercentFormatter(xmax=data['y2'][0], decimals=0))
        plt.savefig("diplomka/pictures/{}_{}.png".format('shell', study))
        plt.show()
if __name__ == '__main__':
    trubka = {'M': 10, 'T1': 45 + 273.15, 'P': 100000, 'Rf': 0.0, 'Medium': [['H2O', 1]]}
    plast = {'M': 2.97, 'T1': 110 + 273.15, 'P': 200000, 'Rf': 0.0, 'Medium': [['H2O', 1]]}
    ostatni = {
        'Uhel': 30.0,
        'Q': 500000.0,
        'MaxL': 7.0,
        'MaxD': 1.2,
        'MaxP': 500000.0,
    }
    param_study = Calculate_param_study(trubka, plast, ostatni)
    # t_p, t_t, h_p, tl, d_out, shell, length
    # vykon, plocha
    '''
    for param in ['t_p', 't_t', 'h_p', 'tl', 'd_out', 'length']:
        param_study.param_study_graph(param, 'vykon')
        param_study.param_study_graph(param, 'plocha')
    param_study.param_study_graph('tl','vykon')
    param_study.param_study_shell_graph('vykon')
    param_study.param_study_shell_graph('plocha')
    '''
    param_study.param_study_graph('tl','vykon')
    param_study.param_study_graph('tl','plocha')