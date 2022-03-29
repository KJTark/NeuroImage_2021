#!/usr/bin/python
# -*- coding: utf-8 -*-

def makeSin_power(p, powpow):
	x = linspace(0, pi, 180)
	pred = sin(x) ** powpow
	pred = pred / max(pred)
	peak = where(pred == 1)[0]

	tmp = abs(x - p['u']).min()
	shiftInd = abs(x - p['u']).argmin()

	pred = roll(pred, shiftInd - peak[0])  # direction of shift is opposite from 'wshift' in MATLAB

	# here - the x axis should be derived based on teh p.cCenters, so that each
	# function peaks properly...then go back to just shifting N times over the
	# space, not 180 times...

	resampInd = ceil(float(len(x)) / float(len(p['x'])))
	pred = pred[0:-1:int(resampInd)]

	return pred

def imagesc(X,scale=None,cmap=None,**kwds):
	""" implements imagesc as in matlab"""

	import pylab

	# pylab.figure()

	if scale is None:
		vmin = None
		vmax = None
	else:
		vmin = scale[0]
		vmax = scale[1]

	if not kwds.has_key('extent'):
		kwds['extent'] = [0, 1, 0, 1]

	if not kwds.has_key('interpolation'):
		kwds['interpolation'] = 'nearest'

	a = pylab.imshow(
		X,
		vmin=vmin,
		vmax=vmax,
		origin='lower',
		cmap=cmap,
		**kwds
		)
	show()
	return a


def load_attributes(attr_file):

	x = os.path.join(attr_file)
	attr = ColumnData(x, header=True)

	 # attr = SampleAttributes(x)

	return attr


def load_nii(nii_file, mask_file, attr):
	"""load experiment dataset"""

	fds = fmri_dataset(samples=os.path.join(nii_file),
						targets=attr.targets, chunks=attr.chunks,
						mask=os.path.join(mask_file))

	return fds


def lag_correction(fds, runTRs, lagTRs):
	"""correct dataset for hemodynamic lag"""

	# split dataset into runs

	nRuns = len(fds) / float(runTRs)
	if int(nRuns) != nRuns:
		print 'Error! number of TRs per run must be a factor of total TRs'
		raise SystemExit

	nRuns = int(nRuns)

	split_fds = []

	for i in range(nRuns):  # split dataset into separate runs
		split_fds.append(fds[i * runTRs:(i + 1) * runTRs])

	# do the shift for each run

	for i in range(len(split_fds)):
		split_fds[i].sa.targets[lagTRs:] = \
		    split_fds[i].sa.targets[:-lagTRs]  # need to shift target labels too

		split_fds[i].sa.censor[lagTRs:] = (split_fds[i]
			.sa.censor[:-lagTRs]) #and censor labels

		#split_fds[i].sa.Acc[lagTRs:] = (split_fds[i]
		#	.sa.Acc[:-lagTRs]) #and censor labels

		#split_fds[i].sa.Match[lagTRs:] = (split_fds[i]
		#	.sa.Match[:-lagTRs])

		split_fds[i] = (split_fds[i])[lagTRs:]

	##  merge back datasets

	fds = split_fds[0]
	for i in range(1, len(split_fds)):
	    fds.append(split_fds[i])

	return fds


if __name__ == '__main__':
	# libraries needed by pymvpa
	import os
	from mvpa2.suite import *
	from numpy import *
	from pylab import *

	nColors = 8  # number of actual orientations in your study
	nChans = 8  # number of channels that you want to use as your 'basis set' to model the response in each voxel.
	# nTrials = 40  # number trials per orientation for the simulation (must be evenly divisible by nColors)

	subs = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','17','19','25']

	for subNo in subs:

		data_dir = subNo + 'Avg_HV_LUM/Img_data_reanal_v2/encoding_AvgHvLum/'
		mask_dir = subNo + 'Avg_HV_LUM/Img_data_reanal_v2/ROI_masks/'
		os.chdir(data_dir) 
		print 'Sub No: %s' % subNo

		attr_file = subNo + 'Avg_HV_LUM/BH_data/' + subNo + '_AvgHvLum_Att_correct_SD7.txt'

		ROIs = ['V1_thr','V2_thr','V3_thr','IPS2_thr','SPL_thr', 'FEF_thr']

		nROIs = 1
		for i in ROIs:

			subjname = (subNo + 'Avg_AvgHvLum_Avg_' + i + '_SD7')

			## ===load files===
			attr = load_attributes(attr_file=attr_file)
			fds = load_nii(nii_file=(subNo + 'Avg_AvgHvLum_Allruns_sc_dt_hp.nii'),
			               mask_file=(mask_dir+i+'.nii'), attr=attr)
			fds.sa['censor'] = attr.censor
			# fds.sa['Acc'] = attr.Acc

			nVox = fds.nfeatures
			print 'Voxel Number Orignial:', nVox
			print 'Initial Number of TRs %s' % len(fds)

			powpow = 6  # sin power for basis function

			fds.samples = asarray(fds.samples)

			## ===lag correction===
			fds = lag_correction(fds=fds, runTRs=76, lagTRs=3)  # another custom subfunction
			print 'after lag corrected: %s' % len(fds)

			## ===censoring===
			fds = fds[fds.sa.censor == 1]  # remove censored points
			print 'after censor removed:%s' % len(fds)

			# ===detrend===
			# poly_detrend(fds, polyord=1, chunks_attr='chunks')

			# ===zscore per run===
			zscore(fds, chunks_attr='chunks')

			# ===remove 'rest' TRs===
			fds = fds[fds.targets != 0]
			print 'after rest TR removed: %s' % len(fds)

			Avg_runs = unique(fds.sa.chunks)
			Avg_runs = Avg_runs[Avg_runs < 99]

			Hv_runs = unique(fds.sa.chunks)
			Hv_runs = Hv_runs[np.logical_and(Hv_runs > 99, Hv_runs < 199)]
			Lum_runs = unique(fds.sa.chunks)
			Lum_runs = Lum_runs[Lum_runs > 199]

			print Avg_runs

			## AVG 
			chan = []
			fds_train_orig = fds[fds.sa.chunks < 99]
			fds_test_orig = fds[fds.sa.chunks < 99]

			#fds_train_orig = fds[fds.sa.chunks < 7]
			#fds_test_orig = fds[fds.sa.chunks < 7]

			for rr in Avg_runs:
				print rr
				
				fds_train = fds_train_orig[fds_train_orig.sa.chunks != rr]
				fds_test = fds_test_orig[fds_test_orig.sa.chunks == rr]

			 	# fds_train = fds_train[fds_train.sa.Acc == 1]
				#fds_test = fds_test[fds_test.sa.Acc == 1]

				# ## feature selection ###
				clf = LinearCSVMC()
				nfeatures = nVox
				
				fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
				        FixedNElementTailSelector(nfeatures, tail='upper',
				        mode='select', sort=False))

				fclf = FeatureSelectionClassifier(clf, fsel)
				fclf.train(fds_train)
				fds_train = fclf.mapper.forward(fds_train)
				fds_test = fclf.mapper.forward(fds_test)
	            # ## end of feature selection ###

				averager = mean_group_sample(['targets', 'chunks'])
				fds_train = fds_train.get_mapped(averager)
				fds_test = fds_test.get_mapped(averager)

				g_train = fds_train.sa.targets
				g_train = g_train.astype(float)
				scan_train = fds_train.sa.chunks
				scan_train = scan_train.astype(float)

				g_test = fds_test.sa.targets
				g_test = g_test.astype(float)
				scan_test = fds_test.sa.chunks
				scan_test = scan_test.astype(float)

				data_train = fds_train.samples  # alloc a matrix for storing data
				temp = data_train[:]
				temp_chunks = fds_train.chunks
				data_test = fds_test.samples  # alloc a matrix for storing data

				data_train = (data_train) + 100
				data_test = (data_test) + 100

				# set up params of the basis function that will govern response of the voxel
				p = {}
				p['x'] = linspace(0, pi - pi / nColors, nColors)  # x-axis to eval the gauss

				print 'Computing iteration {0} out of {1}\n'.format(rr,
				        shape(Hv_runs)[0])

				trn = data_train  # data from training scans (all but one scan)
				tst = data_test  # data from test scan (held out scan)
				trns = scan_train  # vector of scan labels for traning data
				trng = g_train  # vector of trial labels for training data.
				tstg = g_test  # trial labels for tst data.

				uRuns = unique(trns)
				tmp = zeros([nColors * len(uRuns), shape(trn)[1]])
				sCnt = 1
				for ss in uRuns:
					for ii in arange(nColors) + 1:
						tmp[sCnt * nColors - nColors + ii - 1, :] = \
							nanmean(trn[logical_and(trns == ss, trng == ii), :], 0)
					# ....tmp2 = trna[logical_and(trns==ss,trng==ii)]

					sCnt = sCnt + 1
				trn = tmp

				p['x'] = linspace(0, pi - pi / nColors, nColors)
				p['cCenters'] = linspace(0, pi - pi / nChans, nChans)
				X = zeros([shape(trn)[0], nChans])

				for ii in arange(nChans):
					p['u'] = p['cCenters'][ii]
					X[:, ii] = tile(makeSin_power(p, powpow), [1,
									shape(trn)[0] / nColors])

				trn[isnan(trn)] = 0
				X[isnan(X)] = 0
				w = dot(dot(linalg.inv(dot(X.T, X)), X.T), trn)

				x = dot(dot(linalg.inv(dot(w, w.T)), w), tst.T)
				# x = stats.zscore(x)
				x = x.T

				if rr == 1:
					chan = x
					g = g_test
				else:
					chan = vstack([chan, x])
					g = hstack([g, g_test])

			result_each_unshift = zeros([8, 8])
			for k in arange(8):
				result_each_unshift[k, :] = nanmean(chan[g == k + 1], axis=0)

	    	# then shift the rows-data so that the channel corresponding to the stimulus on
	    	# each trial is in the middle column
			for ii in range(shape(chan)[0]):
				chan[ii, :] = roll(chan[ii, :], int(ceil(nChans / 2)
								    - g[ii]))

	         # again, python "roll" is opp direction from matlab "wshift"
			result_each_shift = zeros([8, 8])
			for k in arange(8):
				result_each_shift[k, :] = mean(chan[g == k + 1], axis=0)

			channel_responses_each_unshift = result_each_unshift
			channel_responses_each_shift = result_each_shift


			savetxt(subjname + '_result_all_8ori_unshift.txt', 
				channel_responses_each_unshift, '%f', '\t')
			savetxt(subjname + '_result_all_8ori_shift.txt', 
				channel_responses_each_shift, '%f', '\t')

			## HV ============================================
			subjname = (subNo + 'Avg_AvgHvLum_HV_' + i + '_SD7')
			chan = []
			fds_train = fds[fds.sa.chunks < 99]
			#fds_train = fds[fds.sa.chunks < 7]

			fds_test = fds[np.logical_and(fds.sa.chunks > 99, fds.sa.chunks < 199)]

			for rr in range(0,1):

				# ## feature selection ###
				clf = LinearCSVMC()
				nfeatures = nVox
				
				fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
				        FixedNElementTailSelector(nfeatures, tail='upper',
				        mode='select', sort=False))

				fclf = FeatureSelectionClassifier(clf, fsel)
				fclf.train(fds_train)
				fds_train = fclf.mapper.forward(fds_train)
				fds_test = fclf.mapper.forward(fds_test)
	            # ## end of feature selection ###

				averager = mean_group_sample(['targets', 'chunks'])
				fds_train = fds_train.get_mapped(averager)
				fds_test = fds_test.get_mapped(averager)

				g_train = fds_train.sa.targets
				g_train = g_train.astype(float)
				scan_train = fds_train.sa.chunks
				scan_train = scan_train.astype(float)

				g_test = fds_test.sa.targets
				g_test = g_test.astype(float)
				scan_test = fds_test.sa.chunks
				scan_test = scan_test.astype(float)

				data_train = fds_train.samples  # alloc a matrix for storing data
				temp = data_train[:]
				temp_chunks = fds_train.chunks
				data_test = fds_test.samples  # alloc a matrix for storing data

				data_train = (data_train) + 100
				data_test = (data_test) + 100

				# set up params of the basis function that will govern response of the voxel
				p = {}
				p['x'] = linspace(0, pi - pi / nColors, nColors)  # x-axis to eval the gauss

				print 'Computing iteration {0} out of {1}\n'.format(rr,
				        shape(Hv_runs)[0])

				trn = data_train  # data from training scans (all but one scan)
				tst = data_test  # data from test scan (held out scan)
				trns = scan_train  # vector of scan labels for traning data
				trng = g_train  # vector of trial labels for training data.
				tstg = g_test  # trial labels for tst data.

				uRuns = unique(trns)
				tmp = zeros([nColors * len(uRuns), shape(trn)[1]])
				sCnt = 1
				for ss in uRuns:
					for ii in arange(nColors) + 1:
						tmp[sCnt * nColors - nColors + ii - 1, :] = \
							nanmean(trn[logical_and(trns == ss, trng
									== ii), :], 0)
					# ....tmp2 = trna[logical_and(trns==ss,trng==ii)]

					sCnt = sCnt + 1
				trn = tmp

				p['x'] = linspace(0, pi - pi / nColors, nColors)
				p['cCenters'] = linspace(0, pi - pi / nChans, nChans)
				X = zeros([shape(trn)[0], nChans])

				for ii in arange(nChans):
					p['u'] = p['cCenters'][ii]
					X[:, ii] = tile(makeSin_power(p, powpow), [1,
									shape(trn)[0] / nColors])

				trn[isnan(trn)] = 0
				X[isnan(X)] = 0
				w = dot(dot(linalg.inv(dot(X.T, X)), X.T), trn)

				x = dot(dot(linalg.inv(dot(w, w.T)), w), tst.T)
				# x = stats.zscore(x)
				x = x.T

				if rr == 0:
					chan = x
					g = g_test
				else:
					chan = vstack([chan, x])
					g = hstack([g, g_test])

			result_each_unshift = zeros([8, 8])
			for k in arange(8):
				result_each_unshift[k, :] = nanmean(chan[g == k + 1], axis=0)

	    	# then shift the rows-data so that the channel corresponding to the stimulus on
	    	# each trial is in the middle column
			for ii in range(shape(chan)[0]):
				chan[ii, :] = roll(chan[ii, :], int(ceil(nChans / 2)
								    - g[ii]))

	         # again, python "roll" is opp direction from matlab "wshift"
			result_each_shift = zeros([8, 8])
			for k in arange(8):
				result_each_shift[k, :] = mean(chan[g == k + 1], axis=0)

			channel_responses_each_unshift = result_each_unshift
			channel_responses_each_shift = result_each_shift


			savetxt(subjname + '_result_all_8ori_unshift.txt', 
				channel_responses_each_unshift, '%f', '\t')
			savetxt(subjname + '_result_all_8ori_shift.txt', 
				channel_responses_each_shift, '%f', '\t')



			## Lum  ============================================
			subjname = (subNo + 'Avg_AvgHvLum_Lum_' + i + '_SD7')
			chan = []
			fds_train = fds[fds.sa.chunks < 99]
			fds_test = fds[fds.sa.chunks > 199]

			for rr in range(0,1):
				# fds_test = fds[fds.sa.chunks == rr]
				# fds_train = fds[fds.sa.chunks != rr]
				### example:  fds_cond[np.logical_or(fds_cond.sa.TR == 8, fds_cond.sa.TR == 8+6)]

				
				#fds_train = fds[fds.sa.chunks != Avg_runs[rr]]
				#fds_train = fds_train[fds_train.sa.chunks != Hv_runs[rr]]
				#fds_train = fds_train[fds_train.sa.chunks != Lum_runs[rr]]
				

				# fds_train = fds[np.logical_and(fds.sa.chunks != Avg_runs[rr], fds.sa.chunks != Hv_runs[rr], fds.sa.chunks != Lum_runs[rr])]
				
				# ## feature selection ###
				clf = LinearCSVMC()
				nfeatures = nVox
				
				fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
				        FixedNElementTailSelector(nfeatures, tail='upper',
				        mode='select', sort=False))

				fclf = FeatureSelectionClassifier(clf, fsel)
				fclf.train(fds_train)
				fds_train = fclf.mapper.forward(fds_train)
				fds_test = fclf.mapper.forward(fds_test)
	            # ## end of feature selection ###

				averager = mean_group_sample(['targets', 'chunks'])
				fds_train = fds_train.get_mapped(averager)
				fds_test = fds_test.get_mapped(averager)

				g_train = fds_train.sa.targets
				g_train = g_train.astype(float)
				scan_train = fds_train.sa.chunks
				scan_train = scan_train.astype(float)

				g_test = fds_test.sa.targets
				g_test = g_test.astype(float)
				scan_test = fds_test.sa.chunks
				scan_test = scan_test.astype(float)

				data_train = fds_train.samples  # alloc a matrix for storing data
				temp = data_train[:]
				temp_chunks = fds_train.chunks
				data_test = fds_test.samples  # alloc a matrix for storing data

				data_train = (data_train) + 100
				data_test = (data_test) + 100

				# set up params of the basis function that will govern response of the voxel
				p = {}
				p['x'] = linspace(0, pi - pi / nColors, nColors)  # x-axis to eval the gauss

				print 'Computing iteration {0} out of {1}\n'.format(rr,
				        shape(Lum_runs)[0])

				trn = data_train  # data from training scans (all but one scan)
				tst = data_test  # data from test scan (held out scan)
				trns = scan_train  # vector of scan labels for traning data
				trng = g_train  # vector of trial labels for training data.
				tstg = g_test  # trial labels for tst data.

				uRuns = unique(trns)
				tmp = zeros([nColors * len(uRuns), shape(trn)[1]])
				sCnt = 1
				for ss in uRuns:
					for ii in arange(nColors) + 1:
						tmp[sCnt * nColors - nColors + ii - 1, :] = \
							nanmean(trn[logical_and(trns == ss, trng
									== ii), :], 0)
					# ....tmp2 = trna[logical_and(trns==ss,trng==ii)]

					sCnt = sCnt + 1
				trn = tmp

				p['x'] = linspace(0, pi - pi / nColors, nColors)
				p['cCenters'] = linspace(0, pi - pi / nChans, nChans)
				X = zeros([shape(trn)[0], nChans])

				for ii in arange(nChans):
					p['u'] = p['cCenters'][ii]
					X[:, ii] = tile(makeSin_power(p, powpow), [1,
									shape(trn)[0] / nColors])

				trn[isnan(trn)] = 0
				X[isnan(X)] = 0
				w = dot(dot(linalg.inv(dot(X.T, X)), X.T), trn)

				x = dot(dot(linalg.inv(dot(w, w.T)), w), tst.T)
				# x = stats.zscore(x)
				x = x.T

				if rr == 0:
					chan = x
					g = g_test
				else:
					chan = vstack([chan, x])
					g = hstack([g, g_test])

			result_each_unshift = zeros([8, 8])
			for k in arange(8):
				result_each_unshift[k, :] = nanmean(chan[g == k + 1], axis=0)

	    	# then shift the rows-data so that the channel corresponding to the stimulus on
	    	# each trial is in the middle column
			for ii in range(shape(chan)[0]):
				chan[ii, :] = roll(chan[ii, :], int(ceil(nChans / 2)
								    - g[ii]))

	         # again, python "roll" is opp direction from matlab "wshift"
			result_each_shift = zeros([8, 8])
			for k in arange(8):
				result_each_shift[k, :] = nanmean(chan[g == k + 1], axis=0)

			channel_responses_each_unshift = result_each_unshift
			channel_responses_each_shift = result_each_shift

			savetxt(subjname + '_result_all_8ori_unshift.txt', 
				channel_responses_each_unshift, '%f', '\t')
			savetxt(subjname + '_result_all_8ori_shift.txt', 
				channel_responses_each_shift, '%f', '\t')

			nROIs = nROIs + 1


