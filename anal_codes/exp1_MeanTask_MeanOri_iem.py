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

		split_fds[i].sa.Match[lagTRs:] = (split_fds[i]
			.sa.Match[:-lagTRs])

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

	nColors = 6  # number of actual orientations in your study
	nChans = 6 # number of channels that you want to use as your 'basis set' to model the response in each voxel.

	# nTrials = 24  # number trials per orientation for the simulation (must be evenly divisible by nColors)

	subs = ['01']
	# ROIs = ['V1_thr', 'V2_thr','V3_thr']
	# ROIs = ['IPS_thr', 'SPL_thr', 'FEF_thr']
	# ROIs = ['V1_TarQuad_thr', 'V2_TarQuad_thr','V3_TarQuad_thr']

	for subNo in subs:
		traindata_dir = subNo + '/Img_data/encoding_Tuning/'
		testdata_dir = subNo + '/Img_data/encoding_Mean/'
		mask_dir = subNo + '/Img_data/ROI_masks/'
		os.chdir(testdata_dir) 
		print 'Sub No: %s' % subNo

		train_attr_file = subNo + '_Modelbuilding_All_Att.txt'  
		test_attr_file = subNo + '_Mean_Att.txt'

		nROIs = 1
		nChans= 6
		for i in ROIs:
			subjname = (subNo + 'MeanTask_MeanOri_' + i)

			## ===load train files===
			train_attr = load_attributes(attr_file=(traindata_dir + train_attr_file))
			train_fds = load_nii(nii_file=(traindata_dir + subNo + 'Modelbuilding_sc_dt_hp.nii'),
			               mask_file=(mask_dir + i +'.nii'), attr=train_attr)
			train_fds.sa['censor'] = train_attr.censor
			train_fds.sa['Match'] = train_attr.Match

			## ===load test files===
			test_attr = load_attributes(attr_file=(testdata_dir + test_attr_file))
			test_fds = load_nii(nii_file=(testdata_dir + subNo + 'Main_MeanTask_sc_dt_hp.nii'),
			               mask_file=(mask_dir + i +'.nii'), attr=test_attr)
			test_fds.sa['censor'] = test_attr.censor
			test_fds.sa['Match'] = test_attr.Match


			nVox = train_fds.nfeatures
			print 'Voxel Number Orignial:', nVox
			print 'Initial Number of TRs train %s' % len(train_fds)
			print 'Initial Number of TRs test %s' % len(test_fds)

			powpow = 6  # sin power for basis function

			train_fds.samples = asarray(train_fds.samples)
			test_fds.samples = asarray(test_fds.samples)

			## ===lag correction===
			train_fds = lag_correction(fds=train_fds, runTRs=156, lagTRs=3)  # another custom subfunction
			test_fds = lag_correction(fds=test_fds, runTRs=110, lagTRs=3)  # another custom subfunction
			print 'after lag corrected train: %s' % len(train_fds)
			print 'after lag corrected test: %s' % len(test_fds)

			## ===censoring===
			train_fds = train_fds[train_fds.sa.censor == 1]  # remove censored points
			test_fds = test_fds[test_fds.sa.censor == 1]  # remove censored points
			print 'after censor removed train:%s' % len(train_fds)
			print 'after censor removed test:%s' % len(test_fds)

			# ===zscore per run===
			zscore(train_fds, chunks_attr='chunks')
			zscore(test_fds, chunks_attr='chunks')

			# ===remove 'rest' TRs===
			train_fds = train_fds[train_fds.targets != 0]
			# print 'after rest TR removed train: %s' % len(train_fds)

			test_fds = test_fds[test_fds.targets != 0]
			# print 'after rest TR removed test: %s' % len(test_fds)


			runs = unique(test_fds.sa.chunks)  # find the number of unique runs that we did
			runs = asarray([1])
			chan = []

			for rr in runs:

				fds_train = train_fds
				#fds_test = test_fds[test_fds.sa.chunks == rr]
				fds_test = test_fds

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

	            		## end of feature selection ###
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
				        shape(runs)[0])

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

					sCnt = sCnt + 1
				trn = tmp

				p['x'] = linspace(0, pi - pi / nColors, nColors)
				p['cCenters'] = linspace(0, pi - pi / nChans, nChans)
				X = zeros([shape(trn)[0], nChans])

				for ii in arange(nChans):
					p['u'] = p['cCenters'][ii]
					X[:, ii] = tile(makeSin_power(p, powpow), [1,shape(trn)[0] / nColors])

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

			result_each_unshift = zeros([nChans, nChans])
			for k in arange(nChans):
				result_each_unshift[k, :] = mean(chan[g == k + 1], axis=0)

	    		# then shift the rows-data so that the channel corresponding to the stimulus on
	    		# each trial is in the middle column
			for ii in range(shape(chan)[0]):
				chan[ii, :] = roll(chan[ii, :], int(ceil(nChans / 2)- g[ii]))

	         	# again, python "roll" is opp direction from matlab "wshift"
			result_each_shift = zeros([nChans, nChans])
			for k in arange(nChans):
				result_each_shift[k, :] = mean(chan[g == k + 1], axis=0)

			channel_responses_each_unshift = result_each_unshift
			channel_responses_each_shift = result_each_shift

			savetxt(subjname + '_result_all_unshift.txt', 
				channel_responses_each_unshift, '%f', '\t')
			savetxt(subjname + '_result_all_shift.txt', 
				channel_responses_each_shift, '%f', '\t')

			nROIs = nROIs + 1




