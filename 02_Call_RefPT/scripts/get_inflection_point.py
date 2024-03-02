import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import generic_filter1d
from kneed import KneeLocator
from scipy.optimize import curve_fit

usage = """
Usage:
This script takes a CDT-formatted matrix file from Aggregate Data (single column vector) and determines the inflection point of the histogram of values.

Example: python get_inflection_point.py -i /path/to/MATRIX.cdt -o /path/to/figure.png
"""

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='This script takes a directory of composite (*.out) files and combines them into an excel spreadsheet.')
	parser.add_argument('-i','--input', metavar='composite-dir', required=True, help='directory with all the composite data files')
	parser.add_argument('-o','--output', metavar='outfile', required=True, help='output name to save workbook to')

	args = parser.parse_args()
	return(args)

def loadValues(input):
	'''Load values from CDT (first column)'''
	values = []
	reader = open(input, 'r')
	for line in reader:
		if(line.find("\t")==0):
			continue
		values.append(float( line.strip().split('\t')[1] ))
	reader.close()
	return values

def model(x, a, b, c):
	'''Function to model the curve fit of the tag count frequency bin'''
	return a * np.exp(-b * x) + c

if __name__ == "__main__":
	'''Calculate second derivative of data, plot it against the data, and return first inflection point.'''
	args = getParams()

	sys.stdout.write('Input File: %s\n' % args.input)

	# Load CDT values
	raw_data = loadValues(args.input)

	# Bin the data
	hist_data, bins = np.histogram(raw_data, bins=range(0, int(max(raw_data))))

	# Second derivative approach
	stdv = np.std(raw_data)
	mean = np.mean(raw_data)
	hist_d2 = np.diff(np.diff(hist_data))
	infls = np.where(np.diff(np.sign(hist_d2)))[0]
	sys.stdout.write('Second Derivative: %s\n' % str(infls[0]))

	# # Curve fit and get d2
	# fitopt, fitcov = curve_fit(model, bins[:-1], hist_data)
	# fit_hist = model(bins[:-1], *fitopt)
	# fit_hist_d2 = np.diff(np.diff(fit_hist))
	# fit_infls = np.where(np.diff(np.sign(fit_hist_d2)))[0]
	# sys.stdout.write('Second Derivative (fit): %s\n' % str(fit_infls[0]))

	# KneeLocator approach
	kn = KneeLocator(range(len(hist_data)), hist_data, curve='convex', direction='decreasing')
	sys.stdout.write('KneeLocator: %d\n' % kn.knee)

	sum_total = np.sum(hist_data)
	sum_knee = np.sum(hist_data[kn.knee:])
	sum_d2 = np.sum(hist_data[infls[0]:])
	# sum_fit = np.sum(hist_data[fit_infls[0]:])
	sum_mean = np.count_nonzero(raw_data >= mean)
	sum_stdv = np.count_nonzero(raw_data >= stdv)

	print("Total Sites: %i" % sum_total)
	print("Sites (knee): %i" % sum_knee)
	print("Sites (d2): %i" % sum_d2)
	print("Sites (mean): %i" % sum_mean)
	print("Sites (stdv): %i" % sum_stdv)
	# print("Sites (fit): %i" % sum_fit)

	# Plot raw data, smoothed data, and 1st derivative
	plt.plot(hist_data, label='Raw Data')
	# plt.plot(fit_hist, color='g', label='Fitted Data')

	# Plot inflection points and legend
	plt.axvline(x=kn.knee, color='r', label=f'Inflection Point (knee) {kn.knee} ({sum_knee})')
	plt.axvline(x=stdv, color='b', label=f'Inflection Point (std) {stdv} ({sum_stdv})')
	plt.axvline(x=mean, color='g', label=f'Inflection Point (mean) {mean} ({sum_mean})')
	for i, infl in enumerate(infls, 1):
		sum_infl = np.sum(hist_data[infl:])
		plt.axvline(x=infl, color='k', label=f'Inflection Point (d2) {infl} ({sum_infl})')
		if(i>5):
			break

	plt.xlim(0,50)
	# plt.title("%s - (d2 = %i), (Knee = %i; %i)" % (os.path.basename(args.input).replace("_sorted_SCORES.out",""), infl, kn.knee, sum_knee))
	# plt.legend(bbox_to_anchor=(1,1), bbox_transform=fig.transFigure)
	plt.title(f'{args.input} ({sum_total} sites)')
	plt.legend()

	# Save figure
	plt.savefig(args.output)

	#
	# # Smooth data (check if necessary or even good fit?...maybe necessary to force an inflection...)
	# smooth = gaussian_filter1d(raw_data, 10)
	#
	# # compute second derivative
	# smooth_d2 = np.gradient(np.gradient(smooth))
	#
	# # find switching points
	# infls = np.where(np.diff(np.sign(smooth_d2)))[0]
	#
	# # Plot raw data, smoothed data, and 1st derivative
	# # plt.plot(raw_data, label='Raw Data')
	# # plt.plot(smooth, label='Smoothed Data')
	# # plt.plot(smooth_d2 / np.max(smooth_d2), label='Second Derivative (scaled)')
	# plt.axvline(x=kn.knee, color='k', label=f'Inflection Point {kn.knee}')
	#
	# # Plot inflection points and legend
	# # for i, infl in enumerate(infls, 1):
	# 	# print(infl)
	# 	# plt.axvline(x=infl, color='k', label=f'Inflection Point {i}')
	# plt.legend(bbox_to_anchor=(1.55, 1.0))
	# plt.xlim(0,3000)
	#
	# # Save figure
	# plt.savefig(args.output)
