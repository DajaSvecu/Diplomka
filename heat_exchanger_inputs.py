class UserInputMedium:
    def __init__(self, name):
        self.name = name
        self.parameters = [
            {"name": "M",
                "parse_function": float,
                "default_value": 0.833,
                "unit": "kg/s"},
            {"name": "T1",
                "parse_function": float,
                "default_value": 673.15,
                "unit": "K"},
            {"name": "T2",
                "parse_function": float,
                "default_value": 423.15,
                "unit": "K"},               
            {"name": "P",
                "parse_function": float,
                "default_value": 100000,
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
                "default_value": 100,
                "parse_function": float
            }
        ]

class UserInputShell:
    def __init__(self, name):
        self.name = name
        self.parameters = [
            {"name": "M",
                "parse_function": float,
                "default_value": 10,
                "unit": "kg/s"},
            {"name": "T1",
                "parse_function": float,
                "default_value": 298.15,
                "unit": "K"},
            {"name": "T2",
                "parse_function": float,
                "default_value": 348.15,
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
                "default_value": 100,
                "parse_function": float
            }
        ]

class UserInputRest:
    def __init__(self):
        self.name = 'Rest'
        self.parameters = [
            {
                'name' : 'Uhel',
                "default_value": 30,
                'parse_function': float,
                'unit': 'deg'
            },
            {
                'name' : 'Lambda',
                "default_value": 47,
                'parse_function': float,
                'unit': 'W/m^2K'
            },            
            {
                'name' : 'MaxDelka',
                "default_value": 5,
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxSirka',
                "default_value": 1.2,
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxZtraty',
                "default_value": 50000,
                'parse_function': float,
                'unit': 'Pa'
            }
        ]

class Sizes:
    def __init__(self):
        self.rho = 7850
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
