# this is executed by other scripts that set all the required vars

if 'histogram' in globals() and 'printparams' in globals():
	if histogram and printparams:
		print("Can't choose histogram and printparams in same run.")
		import sys
		sys.exit()
do_select = True
do_process = True
do_play = True
if 'histogram' in globals() and histogram:
	do_select = False
	do_process = False
	do_play = False
	print("Printing histogram of selected partials, before processing. Not playing.")
elif 'printparams' in globals() and printparams:
	do_play = False
	if printparams == 1:
		do_process = False
		print("Printing selected partials, before processing. Not playing.")
	else:
		print("Printing selected partials, after processing. Not playing.")

if (numoutchans == 0 or (numoutchans % 2) != 0):
	print("\nERROR: Must set <numoutchans> to an integer multiple of 2.")
	import sys
	sys.exit()
	
r1 = resynth.Resynth(sr, int(numoutchans / 2))
r2 = resynth.Resynth(sr, int(numoutchans / 2))
print("resynth version:", r1.version())
if 'maxdispargs' in globals():
	r1.set_max_disp_args(maxdispargs)
	r2.set_max_disp_args(maxdispargs)
if 'use_aux_out' in globals() and use_aux_out:
	r1.set_multichan_bustype("aux")
	r2.set_multichan_bustype("aux")

if 'verbosity' in globals() and verbosity > 0:
	if verbosity % 2:
		r1.set_verbose(True)
		r2.set_verbose(True)
	if verbosity > 1:
		rtcmix.set_option("print = 5")

# We always select the freq and time range, even if do_select is false.
#r.select_partial_range(0, 10)	# for debugging
if 'clipend' in globals():
	r1.select_time_range(selstart, selend, clipend)
	r2.select_time_range(selstart, selend, clipend)
else:
	r1.select_time_range(selstart, selend)
	r2.select_time_range(selstart, selend)
if 'minfreq' not in globals():
	minfreq = 20; maxfreq = sr / 2
r1.select_freq_range(minfreq, maxfreq)
r2.select_freq_range(minfreq, maxfreq)

if 'seed' in globals():
	r1.set_pan_seed(seed)
	r2.set_pan_seed(seed)
	r1.set_synth_seed(seed + 1)
	r2.set_synth_seed(seed + 1)
if 'delseed' in globals():
	r1.set_delay_seed(delseed)
	r2.set_delay_seed(delseed)
if 'quantseed' in globals():
	r1.set_quantize_seed(quantseed)
	r2.set_quantize_seed(quantseed)

if do_select:
	if 'mindur' in globals():
		r1.select_duration_range(mindur, maxdur)
		r2.select_duration_range(mindur, maxdur)
	else:
		mindur = 0; maxdur = 10000
	if 'minamp' in globals():
		r1.select_amp_range(minamp, maxamp)
		r2.select_amp_range(minamp, maxamp)
	else:
		minamp = 0; maxamp = 1
	if 'minbw' in globals():
		r1.select_bandwidth_range(minbw, maxbw)
		r2.select_bandwidth_range(minbw, maxbw)
	else:
		minbw = 0; maxbw = 1
	if 'minstd' in globals():
		r1.select_freq_std_range(minstd, maxstd)
		r2.select_freq_std_range(minstd, maxstd)
	else:
		minstd = 0; maxstd = 10000
	if 'fundamental' in globals():
		r1.select_harmonic_partials(fundamental, tolerance)
		r2.select_harmonic_partials(fundamental, tolerance)
	if 'sel_chord' in globals():
		r1.select_chord_partials(sel_chord, sel_transp, sel_sensitivity)
		r2.select_chord_partials(sel_chord, sel_transp, sel_sensitivity)
else:
	# Don't actually select for amp, dur, bw, and std, so that we can get
	# accurate histogram totals for each of these params independently within
	# the time and freq ranges.
	if 'mindur' not in globals():
		mindur = 0; maxdur = 10000
	if 'minamp' not in globals():
		minamp = 0; maxamp = 1
	if 'minbw' not in globals():
		minbw = 0; maxbw = 1
	if 'minstd' not in globals():
		minstd = 0; maxstd = 10000

numpartials1 = r1.read_file(parfile_name1)	# must do this after selects
numpartials2 = r2.read_file(parfile_name2)

if do_process and numpartials1 > 0 and numpartials2 > 0:
	if 'clampminamp' in globals():
		r1.clamp_amps(clampminamp, clampmaxamp)
		r2.clamp_amps(clampminamp, clampmaxamp)
	if 'invert_amps' in globals() and invert_amps:
		r1.invert_amps()
		r2.invert_amps()
	if 'zero_first_amp' in globals() and zero_first_amp:
		r1.set_zero_first_partial_amp(True)
		r2.set_zero_first_partial_amp(True)
	if 'zero_last_amp' in globals() and zero_last_amp:
		r1.set_zero_last_partial_amp(True)
		r2.set_zero_last_partial_amp(True)
	if 'gain' in globals():
		r1.set_gain(gain)
		r2.set_gain(gain)
	if 'clampminbw' in globals():
		r1.clamp_bandwidths(clampminbw, clampmaxbw)
		r2.clamp_bandwidths(clampminbw, clampmaxbw)
	if 'scalebw' in globals() and scalebw != 1.0:
		r1.scale_bandwidths(scalebw)
		r2.scale_bandwidths(scalebw)
	if 'bwcurve' in globals():
		r1.scale_bandwidths(bwcurve)
		r2.scale_bandwidths(bwcurve)
	if 'retune_last' in globals() and retune_last:
		if 'freqscale' in globals():
			r1.scale_and_offset_freqs(freqscale, freqshift, freqdev)
			r2.scale_and_offset_freqs(freqscale, freqshift, freqdev)
		if 'sminfreq' in globals():
			r1.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
			r2.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
		if 'chord' in globals():
			r1.retune_partials(chord, transp, sensitivity, strength)
			r2.retune_partials(chord, transp, sensitivity, strength)
	else:
		if 'chord' in globals():
			r1.retune_partials(chord, transp, sensitivity, strength)
			r2.retune_partials(chord, transp, sensitivity, strength)
		if 'freqscale' in globals():
			r1.scale_and_offset_freqs(freqscale, freqshift, freqdev)
			r2.scale_and_offset_freqs(freqscale, freqshift, freqdev)
		if 'sminfreq' in globals():
			r1.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
			r2.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
	if 'eqtab' in globals():
		r1.eq_partials(eqtab)
		r2.eq_partials(eqtab)
	if 'delaytab' in globals():
		r1.delay_times(delaytab, deldev)
		r2.delay_times(delaytab, deldev)
	if 'quantum' in globals():
		r1.quantize_times(quantum, quantdev, quantsmear, quanttmpl, qfmin, qfmax)
		r2.quantize_times(quantum, quantdev, quantsmear, quanttmpl, qfmin, qfmax)

	if do_play:
		(starttime1, endtime1) = r1.get_time_bounds()
		(starttime2, endtime2) = r2.get_time_bounds()
		starttime = min(starttime1, starttime2)
		endtime = max(endtime1, endtime2)
		print("input time bounds: [{:.6f}, {:.6f}]".format(starttime, endtime))
		timeshift = -starttime
		print("timeshift: {:.6f}".format(timeshift))
		r1.shift_times(timeshift)
		r2.shift_times(timeshift)
		if 'clampdurs' in globals():
			r1.clamp_partial_durations(clampdurs[0], clampdurs[1])
			r2.clamp_partial_durations(clampdurs[0], clampdurs[1])
		if 'timescale_start_points' in globals():
			r1.set_timescale_start_points(timescale_start_points)
			r2.set_timescale_start_points(timescale_start_points)
		if 'additional_timepointscale' in globals() and additional_timepointscale != 1.0:
			r1.scale_time_points(additional_timepointscale)
			r2.scale_time_points(additional_timepointscale)
		if 'timescale' in globals() and timescale != 1.0:
			r1.set_timescale(timescale)
			r2.set_timescale(timescale)
		(starttime1, endtime1) = r1.get_time_bounds()
		(starttime2, endtime2) = r2.get_time_bounds()
		starttime = min(starttime1, starttime2)
		endtime = max(endtime1, endtime2)
		print("output time bounds (pre timescale): [{:.6f}, {:.6f}]".format(starttime, endtime))

		# set other play-time transformations
		if 'envtab' in globals():
			r1.set_ampenv2(envtab)
			r2.set_ampenv2(envtab)
		if 'wavet' in globals():
			r1.set_wavetable(wavet)
			r2.set_wavetable(wavet)
		if 'glisstab' in globals():
			r1.set_glisstable(glisstab, gdur[0], gdur[1])
			r2.set_glisstable(glisstab, gdur[0], gdur[1])
		if 'lfotype' in globals():
			r1.set_lfo(lfotype, lrate, lmin, lmax, lseed, lsmooth, ldmin, ldmax)
			r2.set_lfo(lfotype, lrate, lmin, lmax, lseed, lsmooth, ldmin, ldmax)
		if 'amplfotype' in globals():
			r1.set_amplfo(amplfotype, alrate, almin, almax, alseed, alsmooth, aldmin, aldmax)
			r2.set_amplfo(amplfotype, alrate, almin, almax, alseed, alsmooth, aldmin, aldmax)

		if 'paroutfile_name' in globals():
			print("Writing parfile '{:s}'".format(paroutfile_name1))
			r1.write_file(paroutfile_name1, False)
			print("Writing parfile '{:s}'".format(paroutfile_name2))
			r2.write_file(paroutfile_name2, False)
		else:
			outch = numoutchans
			if 'usereverb' in globals() and usereverb:
				if numoutchans != 2:
					print("\nERROR: If using reverb, must set <numoutchans> to 2.")
					import sys
					sys.exit()
				outch = 2
			if 'writeit' in globals() and writeit:
				rtcmix.set_option("clobber=yes", "play=no")
				rtcmix.rtsetparams(sr, outch)
				rtcmix.rtoutput(outfile_name, "24")
			else:
				rtcmix.rtsetparams(sr, outch)
			rtcmix.control_rate(cntlrate)
			if 'usereverb' in globals() and usereverb:
				rtcmix.load("FREEVERB")
				rtcmix.bus_config("FREEVERB", "aux 0-1 in", "out 0-1")
				r1.set_outbus("aux 0 out")
				(start1, end1) = r1.play_partials()
				r2.set_outbus("aux 1 out")  # must wait until after scheduling left chan
				(start2, end2) = r2.play_partials()
				end = max(end1, end2)
				mix = reverb_mix
				sz = reverb_roomsize
				rd = pow(reverb_roomsize + 1, 2)
				rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
			else:
				if numoutchans == 2:
					r1.set_outbus("out 0")
					r1.play_partials()
					r2.set_outbus("out 1")   # must wait until after scheduling left chan
					r2.play_partials()
				else:
					r1.play_partials()
					r2.play_partials()

if 'histogram' in globals() and histogram:
	r1.print_histogram(minamp, maxamp, mindur, maxdur, minbw, maxbw, minstd, maxstd)
	r2.print_histogram(minamp, maxamp, mindur, maxdur, minbw, maxbw, minstd, maxstd)
elif 'printparams' in globals() and printparams:
	r1.print_freq_range_params(minfreq, maxfreq)
	r2.print_freq_range_params(minfreq, maxfreq)

