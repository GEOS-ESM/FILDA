'''
This module, MCBEF_NAMELIST.py, is designed specifically for handling 
the namelist files used in the Monte Carlo Biphasic Estimation of Fire 
Properties (MCBEF). It provides the necessary tools to parse 
configuration data from text files and structure them into a usable 
format for MCBEF processing.

Main Developer:
    Meng Zhou (MZ)

Supervisor:
    Arlindo da Silva (ADS)

Institutions:
    Global Modeling and Assimilation Office (GMAO),
    the National Aeronautics and Space Administration (NASA)
    Goddard Earth Sciences Technology and Research (GESTAR) II,
    University of Maryland, Baltimore County

Dependencies:
    1	os		Default  Miscellaneous operating system interfaces

'''

import os

def op_list(content):
    return content.replace("'", "").split()

def op_string(content):
    return content

def op_bool(content):
    return True if content.upper() == 'T' else False

def op_int(content):
    return int(content)

def op_float(content):
	# Assuming the content is space-separated, modify if needed
    return float(content)

def op_float_list(content):
	# Assuming the content is space-separated, modify if needed
	return [float(item) for item in content.split()]

def op_int_list(content):
	# Assuming the content is space-separated, modify if needed
	return [int(item) for item in content.split()]

def op_bool_list(content):
	# Assuming the content is space-separated
	return [True if item.upper() == 'T' else False for item in content.split()]

def op_dict(content_lines):
	result = {}
	items = content_lines.split()	
	for item in items:
		if '=' in item:  # Check for the '=' character which separates band and value
			key, value = item.split('=')
			result[key] = float(value)
				
	return result

def parse_keywords_config(file_path):
	keywords = {}
	formatwords = {}
	with open(file_path, 'r') as file:
		for line in file:
			 # Skip empty lines and comments
			if line.strip() and not line.startswith('#'): 
				key, value = line.strip().split(':', 1)
				attr, func = value.split(',')
				func = func.strip()
				if func == "None":
					# Only store the continuation symbol as a direct string
					formatwords[key.strip()] = attr.strip()
				else:
					keywords[key.strip()] = [attr.strip(), eval(func)]
	return keywords, formatwords


# Get the absolute path to the directory where the current 
base_dir = os.path.dirname(os.path.abspath(__file__))
# Append the relative path to the desired directory
config_dir = os.path.join(base_dir, 'config.rc')
# Defined at a module or global level
KEYWORDS, FORMATWORDS = parse_keywords_config(config_dir)

class namelist:
    def __init__(self):
        pass

    def set_attribute(self, key, value):
        setattr(self, key, value)

    def __str__(self):
        attrs = [f"{key}: {getattr(self, key)}" for key in 
                 sorted(self.__dict__.keys())]
        return "\n".join(attrs)

import numpy as np

def namelist_init(input_filename):
	settings = namelist()
	
	print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
	print('   Monte Carlo Biphasic Estimation of Fire Properties (MCBEF) V1.0.1   \n')
	print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
	try:
		with open(input_filename, 'r') as file:
			content_dict = {}
			for line in file:
				print(' ',line.strip())
				for key, (attr, op) in KEYWORDS.items():
				
					if line.startswith(key):
						if key not in content_dict.keys():
							content_dict[key] = ''
							
						content = line.split(':', 1)[1].strip()
						content_dict[key] += ' ' + content
						# Update current key for potential continuations
						current_key = key  
	
				if FORMATWORDS['symb_cont'] in line:
					content = line.split(':', 1)[1].strip()
					content_dict[current_key] += ' ' + content
			
			for key in content_dict.keys():
				content = content_dict[key].strip()
				settings.set_attribute(KEYWORDS[key][0], 
									   KEYWORDS[key][1](content))				
	except FileNotFoundError:
		print(f"Error: The file {input_filename} could not be found.")
	except Exception as e:
		print(f"An error occurred: {e}")
	
	if os.path.isdir(settings.compile_path) == False:
		message = f' - MCBEF NAMELIST: Making {settings.compile_path} '
		print(message)
		os.mkdir(settings.compile_path)
	else:
		message = f' - MCBEF NAMELIST: {settings.compile_path} exists'
		print(message)	
	
	settings.precompile_string=''
	if settings.flag_precompile:
		settings.compile_path = settings.compile_path + 'persistent'
		settings.precompile_string=',compiledir_mode=readonly'
		
# 		if os.path.isdir(settings.compile_path) == False:
# 			message = f' - MCBEF NAMELIST: Making {settings.compile_path} '
# 			print(message)
# 			os.mkdir(settings.compile_path)
# 		else:
# 			message = f' - MCBEF NAMELIST: {settings.compile_path} exists'
# 			print(message)
	
	return settings







