#!/usr/bin/env python3

import numpy as np

from astroquery.ipac.irsa.irsa_dust import IrsaDust
from astroquery.svo_fps import SvoFps

import astropy
import astropy.units as u
from astropy.table import Column, Table
import astropy.coordinates as coord

def get_lambda_4_filter(facility_name, filter_ID, lambda_type='WavelengthRef'):
	filter_list = SvoFps.get_filter_list(facility=facility_name)
	lambda_out = filter_list[filter_list['filterID']==filter_ID][lambda_type]
	return lambda_out

def deredden(src_name, facility_name, filter_ID, val=None,):
	"""
	Function to compute dereddening. If no value provided, will return just A_lambda.
	If val provided, will check units to apply dereddening correctly.

	src_name needs to be searchable by astropy.SkyCoord
	facility_name and filter_ID must match the results found on the webpage:
	http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse
	"""

	src_coord = coord.SkyCoord.from_name(src_name)
	ext_table = IrsaDust.get_extinction_table(src_coord)

	E_BminV_SFD = np.median(ext_table['A_SandF']/ext_table['A_over_E_B_V_SandF'])
	E_BminV_SandF = 0.86*E_BminV_SFD
	R_V = 3.1 #3.07
	A_V_SandF = R_V * E_BminV_SandF

	svo_extinction = Table.read('http://svo2.cab.inta-csic.es/theory/fps/getextlaw.php', 
								format='ascii', units=[u.Angstrom, u.Unit('cm2 / g')], names=['wave','opacity'])
	k_V = 211.4

	lambda_ref = get_lambda_4_filter(facility_name, filter_ID, lambda_type='WavelengthRef')

	A_lambda_over_A_V = np.interp(lambda_ref, svo_extinction['wave'], svo_extinction['opacity'])[0]/k_V

	A_lambda = A_lambda_over_A_V * A_V_SandF

	if val is None:
		return A_lambda

	val *= u.dimensionless_unscaled
	
	if val.unit.is_equivalent(u.Unit('erg / cm / cm / s')):
		print('Assuming your input value is given in erg / (cm2 s)')
		print('If this is incorrect, please provide variable with the correct astropy units')
		print('Your return value is in erg / (cm2 s)')
		return val.to(u.Unit('erg / cm / cm / s')) * 10**(0.4 * A_lambda)
	elif val.unit.is_equivalent(u.Jy):
		print('Assuming your input value is given in Jy')
		print('If this is incorrect, please provide variable with the correct astropy units')
		print('Your return value is in Jy')
		return val.to(u.Jy) * 10**(0.4 * A_lambda)
	elif val.unit.is_equivalent(u.dimensionless_unscaled):
		print('Assuming your input value is given as a magnitude')
		print('If this is incorrect, please provide variable with astropy units')
		print('Your return value is a magnitude')
		return val - A_lambda

if __name__=="__main__":
	src_name = 'B2 1811+31'
	facility_name = 'Sloan'
	filter_ID = 'SLOAN/SDSS.g'
	print(deredden(src_name, facility_name, filter_ID))
