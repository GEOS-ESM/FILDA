import pandas as pd
import os
import numpy as np

class Band:
	def __init__(self, lamda, rsr, transmittance, tau_wvp, 
	             tau_other_gas,band_name):
		self.lamda = lamda
		self.rsr = rsr
		self.transmittance = transmittance
		self.tau_wvp = tau_wvp
		self.tau_other_gas = tau_other_gas
		self.name = band_name


class viirs:
	'''
	Configuration of the VIIRS senor, 
	Read the VIIRS band wavelength, relative sensor response function in 
	µm, and atmospheric transmittance (gas absorption only)
	
	'''
	def __init__(self, rsr_dir: str = None, \
	             band_list: list = ['M11', 'M12', 'M13', 'M14', 
	                                'M15', 'M16', 'DNB', 'I05']):
		"""
		Initializes the viirs object with a directory for sensor 
		responses and a list of bands.
		
		Param:
		rsr_dir: Directory where sensor response files are located.
		band_list: List of band identifiers.
		"""
		if rsr_dir is None:
			# Get the absolute path to the directory where the current 
			# script is located
			base_dir = os.path.dirname(os.path.abspath(__file__))
			# Append the relative path to the desired directory
			rsr_dir = os.path.join(base_dir, 'sensor', 'sensor_viirs')

		self.rsr_dir = rsr_dir
		self.band_list = band_list	
	
		for band in band_list:
			setattr(self, band, {})
		self.read_sensor()
	
	def read_sensor(self):
		for band in self.band_list:
			try:
				file_path = self.rsr_dir + '/' + band + '.csv'
				one_band = pd.read_csv(file_path)
				# Create a Band instance and assign to the corresponding attribute
				band_instance = Band(one_band['lambda'].values, one_band['rsr'].values, 
									 one_band['transmittance'].values, one_band['tau_wvp'].values,
									 one_band['tau_other_gas'].values, band)
				setattr(self, band, band_instance)
			except FileNotFoundError:
				print(f' - Could not find {band} file in {self.rsr_dir}')
				
class GetSensor:
    def __init__(self, v_sensor, sel_bands):
        self.v_sensor = v_sensor
        self.sel_bands = sel_bands
        self.rsr = []
        self.transmittance = []
        self.lambdas = []
        self.tau_wvp = []
        self.tau_other_gas = []

        self.extract_data()

    def extract_data(self):
        for band in self.sel_bands:
            band_info = getattr(self.v_sensor, band, None)
            if band_info is not None:
                self.rsr.append(band_info.rsr)
                self.transmittance.append(band_info.transmittance)
                self.lambdas.append(band_info.lamda)
                self.tau_wvp.append(band_info.tau_wvp)
                self.tau_other_gas.append(band_info.tau_other_gas)
            else:
                print(f"Warning: '{band}' not found in the sensor data. Skipping this band.")

        self.rsr = np.array(self.rsr)
        self.transmittance = np.array(self.transmittance)
        self.lambdas = np.array(self.lambdas)
        self.tau_wvp = np.array(self.tau_wvp)
        self.tau_other_gas = np.array(self.tau_other_gas)

    def __str__(self):
        return (f"Sensor Data:\n"
                f"Lambdas: {self.lambdas}\n"
                f"RSR: {self.rsr}\n"
                f"Transmittance: {self.transmittance}\n"
                f"Tau WVP: {self.tau_wvp}\n"
                f"Tau Other Gas: {self.tau_other_gas}")



				
def get_sensor(v_sensor, sel_bands):
	rsr = []
	transmittance = []
	lambdas = []
	tau_wvp = []
	tau_other_gas = []
	
	for band in sel_bands:
		band_infor = getattr(v_sensor, band, None)
		rsr.append(band_infor.rsr)
		transmittance.append(band_infor.transmittance)
		lambdas.append(band_infor.lamda)
		tau_wvp.append(band_infor.tau_wvp)
		tau_other_gas.append(band_infor.tau_other_gas)
		
	rsr = np.array(rsr)
	lambdas = np.array(lambdas)
	transmittance = np.array(transmittance)
	tau_wvp = np.array(tau_wvp)
	tau_other_gas = np.array(tau_other_gas)    
	
	return lambdas, rsr, transmittance, tau_wvp, tau_other_gas






		
if __name__ == '__main__':
    # Execute when the module is not initialized from an import statement.
	v_sensor = viirs()
	print(v_sensor.DNB.lamda)
	print(dir(v_sensor))				
				
				
				
				
			