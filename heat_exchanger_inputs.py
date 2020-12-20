class UserInputMedium:
    def __init__(self, name):
        self.name = name
        self.parameters = [
            {"name": "M",
                "parse_function": float,
                "unit": "kg/s"},
            {"name": "T1",
                "parse_function": float,
                "unit": "K"},
            {"name": "T2",
                "parse_function": float,
                "unit": "K"},               
            {"name": "P",
                "parse_function": float,
                "unit": "Pa"},               
            {"name": "Rf",
                "parse_function": float,
                "unit": "K*m^2/W"},               
            {"name": "Medium",
                "parse_function": str,
                "unit": "-"}
        ]
class UserInputRest:
    def __init__(self):
        self.name = 'Rest'
        self.parameters = [
            {
                'name' : 'Uhel',
                'parse_function': float,
                'unit': 'deg'
            },
            {
                'name' : 'Lambda',
                'parse_function': float,
                'unit': 'W/m^2K'
            },            
            {
                'name' : 'MaxDelka',
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxSirka',
                'parse_function': float,
                'unit': 'm'
            },
            {
                'name' : 'MaxZtraty',
                'parse_function': float,
                'unit': 'Pa'
            }
        ]