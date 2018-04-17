import csv
import numpy as np
import os
import os.path
import argparse
from matplotlib import pyplot as plt
from scipy.stats import entropy
from enum import Enum

class Distributions(Enum):
	RAYLEIGH = 1
	EXPONENTIAL = 2

parser = argparse.ArgumentParser()
parser.add_argument("--data", type = str )
parser.add_argument("--ch", type = int )
parser.add_argument("--th", type = int )
args = parser.parse_args()

def get_kl_div(data, params,distr):

	rel_size = len(data)
		
	kld = []
	for p in params:
		if(distr == Distributions.RAYLEIGH):
			dist = create_rayleigh(p,rel_size)
		elif(distr == Distributions.EXPONENTIAL):
			dist = create_exponential(p,rel_size)
		temp = entropy(data,dist)
		kld.append(temp)

	min_param_ind = np.argmin(kld)
	min_param = params[min_param_ind]
	print(min_param)
	if(distr == Distributions.RAYLEIGH):
		best_fit = create_rayleigh(min_param,rel_size)
	elif(distr == Distributions.EXPONENTIAL):
		best_fit = create_exponential(min_param,rel_size)

	fig = plt.figure()
	plt.subplot(2,1,1)
	plt.title('Call Duration Distribution')
	plt.xlabel('Call Duration (s)')
	plt.ylabel('Probability density')
	plt.bar(np.arange(len(data)),data,label = "Data")
	plt.plot(best_fit,"r-")

	plt.subplot(2,1,2)
	plt.xlabel('Parameter range')
	plt.ylabel('KL Divergence')
	plt.plot(params,kld)

	fig.savefig('exponential-call-d-150.png')
	plt.show()

def create_exponential(lamda,size):
	dist = []
	t = range(1,size+1)
	for i in t:
		dist.append(lamda*np.exp(-i*lamda))	
	return dist 

def create_rayleigh(sigma,size):
	dist = []
	t = range(1,size+1)
	for i in t:
		dist.append((i/np.square(sigma))*np.exp(-(np.square(i))/(2*np.square(sigma))))	
	return dist

th = int(args.th)
channel_num = int(args.ch)

# Studing each channel [Loop over all channels]
for l in range(1,channel_num+1):
	temp = []
	i = 0
	flag = 0
	bin_width = 5
	count = 0
	offset = 5

	occ_time = 0
	call_d = []
	t_start = 0
	with open(args.data.strip(),'rb') as csvfile:
		data = csv.reader(csvfile, delimiter=',')
		for row in data:
			if(float(row[offset + l]) > th and flag == 0):
				temp.append(i)	# Flag the call start time
				t_start = i
				flag = 1
			elif(float(row[offset + l]) > th):
				occ_time+=1
			if(float(row[offset + l]) < th):
				flag = 0
				t_end = i
				call_d.append(t_end - t_start)
				# Flag the call end time.
				# temp2.append(i)
			i = i + 1

	totaltime = i
	print("Number of calls for channel " + str(l) + ": " + str(len(temp)))
	print("Occ Time for channel " + str(l) + ": " + str(occ_time) + " seconds or "+ str(occ_time/60) + " mins.")
	print("Total time: " + str(i) + " seconds.")
	pocc = (occ_time*100 / totaltime)
	print("Percentage occupancy for channel " + str(l) + ": " + str(pocc)+ "%")
	# Calculating Inter-Arrival: time between two consecutive calls 
	iat = [temp[0]]
	for i in range(1,len(temp)):
		iat.append(temp[i] - temp[i-1])

	plt_bin_width = 1
	
	# Cutting down tail of distribution
	rel_size = 40 # for IAT
	# rel_size = 150 # for call_duration
	
	# bins_call = np.arange(0,len(call_d)/plt_bin_width,plt_bin_width)
	bins = np.arange(0, rel_size, plt_bin_width) # Stripping the observed data to relevent size

	# Call Duration histogram creation
	call_d = np.array(call_d)
	call_d_hist,bin_edges_call = np.histogram(call_d,bins)
	call_d_hist_norm = (call_d_hist/float(sum(call_d_hist))).astype(float)

	# Pre-processing data for "0" bins
	for i in range(len(call_d_hist_norm)):
		if(call_d_hist_norm[i] == 0):
			call_d_hist_norm[i] = 0.000001

	# IAT histogram creation
	iat = np.array(iat)
	iat_hist,bin_edges = np.histogram(iat,bins)
	iat_hist_norm = (iat_hist/float(sum(iat_hist))).astype(float)

	# Pre-processing data for "0" bins
	for i in range(len(iat_hist_norm)):
		if(iat_hist_norm[i] == 0):
			iat_hist_norm[i] = 0.000001

	sigma = np.arange(8,20,0.01)
	lamda = np.arange(0.01,0.1,0.001)

	get_kl_div(iat_hist_norm,lamda,Distributions.EXPONENTIAL)

	# get_kl_div(iat_hist_norm,sigma,Distributions.RAYLEIGH)

	# get_kl_div(call_d_hist_norm,lamda,Distributions.EXPONENTIAL)

	# get_kl_div(call_d_hist_norm,sigma,Distributions.RAYLEIGH)
