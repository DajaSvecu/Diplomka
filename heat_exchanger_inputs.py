class UserInputMedium:
    def __init__(self, name, title):
        self.name = name
        self.title = title
        self.parameters = [
            {"name": "M",
                "parse_function": float,
                "default_value": 720,
                "to_SI" : 1/3600,
                "hint" : "Prutok",
                "unit": "kg/hod"},
            {"name": "T1",
                "parse_function": float,
                "default_value": 453.15,
                "to_SI" : 1,
                "hint" : "Vstupni teplota",
                "unit": "K"},
            {"name": "T2",
                "parse_function": float,
                "default_value": 0,
                "to_SI" : 1,
                "hint" : "Vystupni teplota",
                "unit": "K"},               
            {"name": "P",
                "parse_function": float,
                "default_value": 100,
                "to_SI" : 1000,
                "hint" : "Tlak media",
                "unit": "kPa"},               
            {"name": "Rf",
                "parse_function": float,
                "default_value": 0.0001,
                "to_SI" : 1,
                "hint" : "Zanaseni",
                "unit": "K*m^2/W"}
        ]
        self.medium =[
            {
                "name": 'MEDIUM',
                "default_value": 'H2O',
                "to_SI" : 1,
                "hint" : "Chemicka znacka media (H2O, CO2...)",
                "parse_function": str
            },
            {
                "name": 'PROCENTA',
                "default_value": 100,
                "to_SI" : 0.01,
                "hint" : "Procentualni zastoupeni",
                "parse_function": float
            }
        ]

class UserInputShell:
    def __init__(self, name):
        self.name = name
        self.parameters = [
            {"name": "M",
                "parse_function": float,
                "default_value": 1,
                "unit": "kg/hod"},
            {"name": "T1",
                "parse_function": float,
                "default_value": 298.15,
                "unit": "K"},
            {"name": "T2",
                "parse_function": float,
                "default_value": 0,
                "unit": "K"},               
            {"name": "P",
                "parse_function": float,
                "default_value": 400000,
                "unit": "Pa"},               
            {"name": "Rf",
                "parse_function": float,
                "default_value": 0.0001,
                "unit": "K*m^2/W"}
        ]
        self.medium =[
            {
                "name": 'NAZEV',
                "default_value": 'H2O',
                "parse_function": str
            },
            {
                "name": 'PERCENT',
                "default_value": 1,
                "parse_function": float
            }
        ]

class UserInputRest:
    def __init__(self):
        self.name = 'Rest'
        self.title = 'Ostatni'
        self.parameters = [
            {
                'name' : 'Q',
                "default_value": 20,
                "to_SI" : 1000,
                "hint" : "Vykon vymeniku",
                'parse_function': float,
                'unit': 'kW'
            },
            {
                'name' : 'Uhel',
                "default_value": 30,
                "to_SI" : 1,
                "hint" : "Uhel mezi trubkami",
                'parse_function': float,
                'unit': 'deg'
            },            
            {
                'name' : 'MaxL',
                "default_value": 5,
                "to_SI" : 1,
                "hint" : "Maximalni delka vymeniku",
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxD',
                "default_value": 1.2,
                "to_SI" : 1,
                "hint" : "Maximalni prumer plaste vymeniku",
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxP',
                "default_value": 50,
                "to_SI" : 1000,
                "hint" : "Maximalni tlakove ztraty vymeniku",
                'parse_function': float,
                'unit': 'kPa'
            }
        ]

class Sizes:
    def __init__(self):
        self.rho = 7850
        self.lamb = 50
        '''
        self.list_of_tubes = [
            {
                'd_in' : 0.012,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003],
            },
            {
                'd_in' : 0.014,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003],
            },
            {
                'd_in' : 0.016,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003],
            },
            {
                'd_in' : 0.018,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003],
            },
            {
                'd_in' : 0.02,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004],
            },
            {
                'd_in' : 0.022,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004],
            },
            {
                'd_in' : 0.024,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004],
            },
            {
                'd_in' : 0.026,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004],
            },
            {
                'd_in' : 0.028,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004],
            },
            {
                'd_in' : 0.03,
                'wall' : [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005],
            },
        ]'''
        
        self.list_of_tubes = [
            {
                'd_out' : 0.00635,
                'wall' : [0.000889, 0.001244, 0.00165],
            },
            {
                'd_out' : 0.009525,
                'wall' : [0.000889, 0.001244, 0.00165, 0.002108],
            },
            {
                'd_out' : 0.0127,
                'wall' : [0.000889, 0.001244, 0.00165, 0.002108],
            },
            {
                'd_out' : 0.015875,
                'wall' : [0.001244, 0.00165, 0.002108, 0.00277],
            },
            {
                'd_out' : 0.01905,
                'wall' : [0.00165, 0.002108, 0.00277],
            },
            {
                'd_out' : 0.022225,
                'wall' : [0.00165, 0.002108, 0.00277],
            },
            {
                'd_out' : 0.0254,
                'wall' : [0.00165, 0.002108, 0.00277],
            },
            {
                'd_out' : 0.03175,
                'wall' : [0.00165, 0.002108, 0.00277],
            },
            
        ]


        self.list_of_shells = [
            {
                'DN': 200, # [-]
                'wall': 0.008179, # [m]
                'D2': 0.219, # [m]
            },
            {
                'DN': 250, # [-]
                'wall': 0.009271, # [m]
                'D2': 0.273, # [m]
            },
            {
                'DN': 300, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.3239, # [m]
            },
            {
                'DN': 350, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.3556, # [m]
            },
            {
                'DN': 400, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.4064, # [m]
            },
            {
                'DN': 450, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.4572, # [m]
            },
            {
                'DN': 500, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.508, # [m]
            },
            {
                'DN': 550, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.5588, # [m]
            },
            {
                'DN': 600, # [-]
                'wall': 0.009525, # [m]
                'D2': 0.6096, # [m]
            },
            {
                'DN': 650, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.6604, # [m]
            },
            {
                'DN': 700, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.7112, # [m]
            },
            {
                'DN': 750, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.762, # [m]
            },
            {
                'DN': 800, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.8128, # [m]
            },
            {
                'DN': 850, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.8636, # [m]
            },
            {
                'DN': 900, # [-]
                'wall': 0.0127, # [m]
                'D2': 0.9144, # [m]
            },
            {
                'DN': 1000, # [-]
                'wall': 0.0127, # [m]
                'D2': 1.016, # [m]
            },            
        ]
        
        self.list_of_lengths = {1.2192, 1.829, 2.438, 3.048, 3.658, 4.267, 4.877, 5.486, 6.096}

        self.list_of_pitches = {1.25, 1.33, 1.5}
        
        self.list_of_baffle_spaces = {0.3, 0.6, 0.9}