#!/usr/bin/env python3

import numpy as np

import configparser

from astroquery.ipac.irsa.irsa_dust import IrsaDust
from astroquery.svo_fps import SvoFps

import astropy
import astropy.units as u
from astropy.table import Column, Table, QTable
import astropy.coordinates as coord
u.set_enabled_equivalencies(u.spectral())

def get_lambda_4_filter(facility_name, filter_ID, lambda_type='WavelengthRef'):
	filter_list = SvoFps.get_filter_list(facility=facility_name)
	lambda_out = filter_list[filter_list['filterID']==filter_ID][lambda_type].quantity[0]
	return lambda_out

def deredden(src_name, facility_name, filter_ID, val=None, convertJy2Eflux=False):
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

	svo_extinction = QTable.read('http://svo2.cab.inta-csic.es/theory/fps/getextlaw.php', 
								format='ascii', units=[u.Angstrom, u.Unit('cm2 / g')], names=['wave','opacity'])
	k_V = 211.4 * u.Unit('cm2 / g')

	lambda_ref = get_lambda_4_filter(facility_name, filter_ID, lambda_type='WavelengthRef')

	A_lambda_over_A_V = np.interp(lambda_ref, svo_extinction['wave'], svo_extinction['opacity'])/k_V

	A_lambda = A_lambda_over_A_V * A_V_SandF

	if val is None:
		return A_lambda

	val *= u.dimensionless_unscaled
	old_unit = val.unit	
	if val.unit.is_equivalent(u.Unit('erg / (cm2 s)')):
		print('Assuming your input value is given in erg / (cm2 s)')
		print('If this is incorrect, please provide variable with the correct astropy units')
		print('Your return value is in erg / (cm2 s)')
		return val.to(u.Unit('erg / (cm2 s)')) * 10**(0.4 * A_lambda)
	elif val.unit.is_equivalent(u.Jy):
		print('Assuming your input value is given as a spectral flux density')
		print('If this is incorrect, please provide variable with the correct astropy units')
		if convertJy2Eflux:
			print('Your return value is in erg / (cm2 s)')
			return (val.to(u.Jy) * lambda_ref.to(u.Hz) * 10**(0.4 * A_lambda)).to(u.Unit('erg / (cm2 s)'))
		else:
			print('Your return value is given as a spectral flux density')
			return (val.to(u.Jy) * 10**(0.4 * A_lambda)).to(old_unit)
	elif val.unit.is_equivalent(u.dimensionless_unscaled):
		print('Assuming your input value is given as a magnitude')
		print('If this is incorrect, please provide variable with astropy units')
		print('Your return value is a magnitude')
		return val - A_lambda

if __name__=="__main__":
	config = configparser.ConfigParser()
	config.read('config_deredden.ini')

	src_name = config['src']['src_name']
	print(f'Computing dereddening for source: {src_name}')


	sed_in_fn = config['io']['data_file']
	dat = QTable.read(sed_in_fn, format='csv')
	dat['F(mJy)'].unit = u.mJy
	dat['F_errn(mJy)'].unit = u.mJy
	dat['F_errp(mJy)'].unit = u.mJy
	print(dat)

	out_table = QTable(names=['nu', 'nu_min', 'nu_max', 'flux', 'flux_errn', 'flux_errp', 'filter'],
						units=[u.Hz, u.Hz, u.Hz, 
								u.Unit('erg / (cm2 s)'), u.Unit('erg / (cm2 s)'), u.Unit('erg / (cm2 s)'),None],
						dtype=[None]*6+['S2'])

	for row in dat:
		flux_extcorr, flux_errn_extcorr, flux_errp_extcorr = deredden(src_name, row['facility_name'], row['filter_ID'], 
																		val=u.Quantity([row['F(mJy)'], row['F_errn(mJy)'], row['F_errp(mJy)']]), 
																		convertJy2Eflux=True)
		#flux_extcorr = deredden(src_name, row['facility_name'], row['filter_ID'], val=row['F(mJy)'], convertJy2Eflux=True)
		#flux_errn_extcorr = deredden(src_name, row['facility_name'], row['filter_ID'], val=row['F_errn(mJy)'], convertJy2Eflux=True)
		#flux_errp_extcorr = deredden(src_name, row['facility_name'], row['filter_ID'], val=row['F_errp(mJy)'], convertJy2Eflux=True)
		lambda_ref = get_lambda_4_filter(row['facility_name'], row['filter_ID'], lambda_type='WavelengthRef')
		nu = lambda_ref.to(u.Hz)
		fwhm = get_lambda_4_filter(row['facility_name'], row['filter_ID'], lambda_type='FWHM')
		lambda_ref_err = fwhm / np.sqrt(8 * np.log(2))
		nu_err = ((lambda_ref-lambda_ref_err).to(u.Hz) - (lambda_ref+lambda_ref_err).to(u.Hz)) / 2. 
		#nu_err = lambda_ref_err.to(u.Hz)
		out_table.add_row([nu, nu_err, nu_err, flux_extcorr, flux_errn_extcorr, flux_errp_extcorr, row['filter_ID']])
	
	for col in out_table.itercols():
		if col.info.dtype.kind == 'f':
			col.info.format = '{:.4e}'

	print(out_table)
	out_table.write(f"results/{config['io']['output_file']}", format='ascii.ecsv', overwrite=True)

	#facility_name = 'Sloan'
	#filter_ID = 'SLOAN/SDSS.g'
	#print(deredden(src_name, facility_name, filter_ID))
