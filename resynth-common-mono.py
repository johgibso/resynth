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

r = resynth.Resynth(sr, numoutchans)
print("resynth version:", r.version())
if 'maxdispargs' in globals():
	r.set_max_disp_args(maxdispargs)
if 'use_aux_out' in globals() and use_aux_out:
	r.set_multichan_bustype("aux")

if 'verbosity' in globals() and verbosity > 0:
	if verbosity % 2:
		r.set_verbose(True)
	if verbosity > 1:
		rtcmix.set_option("print = 5")

# We always select the freq and time range, even if do_select is false.
#r.select_partial_range(0, 10)	# for debugging
if 'clipend' in globals():
	r.select_time_range(selstart, selend, clipend)
else:
	r.select_time_range(selstart, selend)
if 'minfreq' not in globals():
	minfreq = 20; maxfreq = sr / 2
r.select_freq_range(minfreq, maxfreq)

if 'seed' in globals():
	r.set_pan_seed(seed)
	r.set_synth_seed(seed + 1)
if 'delseed' in globals():
	r.set_delay_seed(delseed)
if 'quantseed' in globals():
	r.set_quantize_seed(quantseed)

if do_select:
	if 'mindur' in globals():
		r.select_duration_range(mindur, maxdur)
	else:
		mindur = 0; maxdur = 10000
	if 'minamp' in globals():
		r.select_amp_range(minamp, maxamp)
	else:
		minamp = 0; maxamp = 1
	if 'minbw' in globals():
		r.select_bandwidth_range(minbw, maxbw)
	else:
		minbw = 0; maxbw = 1
	if 'minstd' in globals():
		r.select_freq_std_range(minstd, maxstd)
	else:
		minstd = 0; maxstd = 10000
	if 'fundamental' in globals():
		r.select_harmonic_partials(fundamental, tolerance)
	if 'sel_chord' in globals():
		r.select_chord_partials(sel_chord, sel_transp, sel_sensitivity)
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

numpartials = r.read_file(parfile_name)	# must do this after selects

if do_process and numpartials > 0:
	if 'clampminamp' in globals():
		r.clamp_amps(clampminamp, clampmaxamp)
	if 'invert_amps' in globals() and invert_amps:
		r.invert_amps()
	if 'zero_first_amp' in globals() and zero_first_amp:
		r.set_zero_first_partial_amp(True)
	if 'zero_last_amp' in globals() and zero_last_amp:
		r.set_zero_last_partial_amp(True)
	if 'gain' in globals():
		r.set_gain(gain)
	if 'clampminbw' in globals():
		r.clamp_bandwidths(clampminbw, clampmaxbw)
	if 'scalebw' in globals() and scalebw != 1.0:
		r.scale_bandwidths(scalebw)
	if 'bwcurve' in globals():
		r.scale_bandwidths(bwcurve)
	if 'retune_last' in globals() and retune_last:
		if 'freqscale' in globals():
			r.scale_and_offset_freqs(freqscale, freqshift, freqdev)
		if 'sminfreq' in globals():
			r.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
		if 'chord' in globals():
			r.retune_partials(chord, transp, sensitivity, strength)
	else:
		if 'chord' in globals():
			r.retune_partials(chord, transp, sensitivity, strength)
		if 'freqscale' in globals():
			r.scale_and_offset_freqs(freqscale, freqshift, freqdev)
		if 'sminfreq' in globals():
			r.synthesize_new_partials(sminfreq, smaxfreq, smult, sgain, sdev)
	if 'eqtab' in globals():
		r.eq_partials(eqtab)
	if 'delaytab' in globals():
		r.delay_times(delaytab, deldev)
	if 'quantum' in globals():
		r.quantize_times(quantum, quantdev, quantsmear, quanttmpl, qfmin, qfmax)

	if do_play:
		(starttime, endtime) = r.get_time_bounds()
		print("input time bounds: [{:.6f}, {:.6f}]".format(starttime, endtime))
		timeshift = -starttime
		print("timeshift: {:.6f}".format(timeshift))
		r.shift_times(timeshift)
		if 'clampdurs' in globals():
			r.clamp_partial_durations(clampdurs[0], clampdurs[1])
		if 'timescale_start_points' in globals():
			r.set_timescale_start_points(timescale_start_points)
		if 'additional_timepointscale' in globals() and additional_timepointscale != 1.0:
			r.scale_time_points(additional_timepointscale)
		if 'timescale' in globals() and timescale != 1.0:
			r.set_timescale(timescale)
		(starttime, endtime) = r.get_time_bounds()
		print("output time bounds (pre timescale): [{:.6f}, {:.6f}]".format(starttime, endtime))

		# set other play-time transformations
		if 'envtab' in globals():
			r.set_ampenv2(envtab)
		if 'wavet' in globals():
			r.set_wavetable(wavet)
		if 'glisstab' in globals():
			r.set_glisstable(glisstab, gdur[0], gdur[1])
		if 'lfotype' in globals():
			r.set_lfo(lfotype, lrate, lmin, lmax, lseed, lsmooth, ldmin, ldmax)
		if 'amplfotype' in globals():
			r.set_amplfo(amplfotype, alrate, almin, almax, alseed, alsmooth, aldmin, aldmax)

		if 'paroutfile_name' in globals():
			print("Writing parfile '{:s}'".format(paroutfile_name))
			r.write_file(paroutfile_name, False)
		else:
			if 'usereverb' in globals() and usereverb:
				if (numoutchans != 1 and numoutchans != 2 and numoutchans != 4 and numoutchans != 8):
					print("\nERROR: If using reverb, must set <numoutchans> to 1, 2, 4, or 8.")
					import sys
					sys.exit()
			if 'writeit' in globals() and writeit:
				rtcmix.set_option("clobber=yes", "play=no")
				rtcmix.rtsetparams(sr, numoutchans)
				rtcmix.rtoutput(outfile_name, "24")
			else:
				rtcmix.rtsetparams(sr, numoutchans)
			rtcmix.control_rate(cntlrate)
			if 'usereverb' in globals() and usereverb:
				rtcmix.load("FREEVERB")
				mix = reverb_mix
				sz = reverb_roomsize
				rd = pow(reverb_roomsize + 1, 2)
				if numoutchans == 1:
					r.set_outbus("aux 0 out")
					(start, end) = r.play_partials()
					rtcmix.bus_config("FREEVERB", "aux 0 in", "out 0-1")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
				elif numoutchans == 2:
					r.set_outbus("aux 0-1 out")
					(start, end) = r.play_partials()
					rtcmix.bus_config("FREEVERB", "aux 0-1 in", "out 0-1")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
				elif numoutchans == 4:
					r.set_multichan_bustype("aux")
					(start, end) = r.play_partials()
					rtcmix.bus_config("FREEVERB", "aux 0-1 in", "out 0-1")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
					rtcmix.bus_config("FREEVERB", "aux 2-3 in", "out 2-3")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
				elif numoutchans == 8:
					r.set_multichan_bustype("aux")
					(start, end) = r.play_partials()
					rtcmix.bus_config("FREEVERB", "aux 0-1 in", "out 0-1")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
					rtcmix.bus_config("FREEVERB", "aux 2-3 in", "out 2-3")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
					rtcmix.bus_config("FREEVERB", "aux 4-5 in", "out 4-5")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
					rtcmix.bus_config("FREEVERB", "aux 6-7 in", "out 6-7")
					rtcmix.FREEVERB(0, 0, end, 1, sz, 0.1, rd, 0, 100-mix, mix, 100)
			else:
				if numoutchans == 1:
					r.set_outbus("out 0")
				r.play_partials()

if 'histogram' in globals() and histogram:
	r.print_histogram(minamp, maxamp, mindur, maxdur, minbw, maxbw, minstd, maxstd)
elif 'printparams' in globals() and printparams:
	r.print_freq_range_params(minfreq, maxfreq)

