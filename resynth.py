# Read a sinusoidal analysis file written by SPEAR in par-text-partials-format
# (File > Export Format > Text - Partials), and perform oscillator resynthesis
# of the file, with altered params, using the RTcmix Python front end, PYCMIX,
# and the RTcmix instrument, BWESINE. Will correctly resynthesize analyses
# done with "time reassignment transient sharpening" option turned on when
# analyzing in SPEAR.
#
# This module will also read files exported from Loris and converted via
# a Python module, "lorispartext.py" (included here), that writes
# par-text-partials-format and a variant of this format that also includes the
# bandwidth data, partial starting phase, and partial label gathered from a
# Loris analysis. 
#
# Make this resynth module visible to Python somehow (by simply writing a
# sym link to it into the directory that contains your scripts), then run a
# script that uses the Resynth class as follows:
#
# pycmix < myscript.py
#
# See "resynth-tmpl-mono-in.py" for an example. This requires the file
# "resynth-common-mono.py," which includes helper code that is likely to be
# used commonly by your scripts and which must reside in the same folder as
# your script (unless you change the path to it in your script). This script
# can write stereo or multichannel output. There are also stereo-input versions
# of these scripts.
#
# NOTE: It's best to recompile RTcmix with MAXDISPARGS set to 4096 (in
# rtcmix/src/include/maxdispargs.h). If you don't do this, files containing
# partials with more than about 500 breakpoints will not work. MAXDISPARGS
# should be a power of 2 a little greater than twice the maximum number of 
# breakpoints you want to use. You might also consider increasing MAXCHANS
# to something greater than the default 8 (in rtcmix/site.conf), if you want
# to work with large spatial audio systems.
#
# John Gibson, 1/15/19
#
#
# =============================================================================
# public class methods
#
# (See the source code for more complete usage comments.)
#
# First, use any of these methods to filter out partials read from file.
#
#    select_partial_range(min, max=0)
#    select_time_range(min, max=0, clipend=False)
#    select_duration_range(min, max=0)
#    select_amp_range(min, max=0)
#    select_freq_range(min, max=0)
#    select_freq_std_range(min, max=0))
#    select_harmonic_partials(fundamental, tolerance=0)
#    select_chord_partials(chord, transp=0, sensitivity=6)
#    select_bandwidth_range(min=0, max=1)
#
# Then read the file.
#
#    numpartials = read_file(filename)
#
# If numpartials > 0, alter partials with the following methods. Ignore the
# partials argument and the return value to affect all the filtered partials
# read by the last call to read_file. These methods can be called in any order,
# possibly with different results.
#
#    partials = retune_partials(chord, transp=0, sensitivity=6, strength=0.9, partials=None)
#    partials = scale_and_offset_freqs(scale, offset=0, freqdev=0, fmin=0, fmax=0, partials=None)
#    partials = eq_partials(eqtable, fmin=0, fmax=0, partials=None)
#    partials = clamp_amps(min, max, partials=None)
#    partials = invert_amps(partials=None)
#    partials = clamp_bandwidths(min, max, partials=None)
#    partials = scale_bandwidths(scale, partials=None)
#    partials = shift_times(shift, partials=None)
#    partials = delay_times(delay, deldev=0, fmin=0, fmax=0, partials=None)
#    partials = quantize_times(quantum, quantdev=0, smear=0, template=0, fmin=0, fmax=0, partials=None)
#    partials = scale_time_points(scale, partials=None)
#    partials = synthesize_new_partials(srcminfreq, srcmaxfreq, freqscale, gain=0, freqdev=0, minfreq=20, partials=None)
#
# Before playing partials, optionally set up a few things.
#
#    set_outbus(bus)
#    set_multichan_bustype(bustype)
#    set_gain(gain)
#    set_zero_first_partial_amp(doit)
#    set_zero_last_partial_amp(doit)
#    set_timescale(timescale)
#    set_timescale_start_points(doit)
#    set_wavetable(table_handle)
#    set_ampenv2(table_handle)
#    set_glisstable(table_handle, mindur=0, maxdur=0)
#    set_lfo(type, rate, min, max, seed=1, smooth=0, mindur=0, maxdur=0)
#    set_amplfo(type, rate, min, max, seed=1, smooth=0, mindur=0, maxdur=0)
#    set_zero_phase(doit)
#    set_pan_seed(seed)
#    set_synth_seed(seed)
#    set_delay_seed(seed)
#    set_quantize_seed(seed)
#    set_gliss_seed(seed)
#    set_lfo_seed(seed)
#    set_amplfo_seed(seed)
#    set_max_disp_args(max)
#
# Then play the partials.
#
#    (minstart, maxend) = play_partials(partials=None)
#
# Some additional methods...
#
#    version_string = version()
#    partials = get_partials()
#    (starttime, endtime) = get_time_bounds(partials=None)
#    (minfreq, maxfreq) = get_partial_freq_range(partial_breakpoints)
#    mean = get_partial_freq_mean(partial_breakpoints, weighted=True)
#    std = get_partial_freq_std(partial_breakpoints, mean=None)
#    (minamp, maxamp) = get_partial_amp_range(partial_breakpoints)
#    (minbw, maxbw) = get_partial_bw_range(partial_breakpoints)
#    set_verbose(verbose)
#    result = write_file(filename, extended=False, partials=None)
#    print_histogram(minamp, maxamp, mindur, maxdur, minbw, maxbw, minstd, maxstd, partials=None)
#    print_freq_range_params(minfreq, maxfreq=0, partials=None)
#
# -----------------------------------------------------------------------------
# Non-class utility functions
#
#    pitchlist = pclist2pitchlist(pclist, minoctave=4, maxoctave=14)
#
#
# =============================================================================
# TODO
# - filters affecting frequency scaling. IOW, scale only frequencies that fit
#   certain criteria. Maybe best to have a general framework that could be
#   used for time-scaling, freq-shifting also? Not sure partial groups approach
#   in my resynth~ Max external is fruitful enough.
#

import copy, math, random, rtcmix, sys, types

class Resynth:
	_kEnvSlotsPerPoint = 2

	def __init__(self, sampling_rate, num_chans):
		self.set_max_disp_args()
		self._version = "0.0.1"
		self._sr = sampling_rate
		self._nyquist = sampling_rate / 2.0
		self._nchans = num_chans
		self._multichan_bustype = "out"
		self._verbose = False
		self._extended = False	# reading our extended partext format
		self._partial_range = (0, 0)
		self._time_range = (0, 0)
		self._clip_end = False
		self._duration_range = (0, 0)
		self._amp_range = (0, 1)
		self._freq_range = (20, self._nyquist)
		self._freq_std_range = (0, 10000)
		self._freq_harmonics = None
		self._freq_harm_tolerance = 0
		self._freq_chord_pitches = None
		self._freq_chord_transp = 0
		self._freq_chord_sensitivity = 0
		self._bw_range = (0, 1)
		self._totdur = 0		# total duration spanned by unfiltered partials
		self._selected_partials_start = 0 # earliest start time, not necc. 1st partial
		self._ampscale = 1.0
		self._retune_chord = []
		self._retune_chord_transp = 0
		self._numretunedpartials = 0
		self._retune_sensitivity = 6.0
		self._retune_strength = 0.9
		self._timescale = 1
		self._timescale_start_points = True
		self._wavetable = rtcmix.maketable("wave", 32767, "sine")
		self._ampenv2 = None
		self._glisstable = None
		self._glissdurrange = (0, 0)
		self._lfotype = None
		self._lfotype_israndom = False
		self._lforate = 5.0
		self._lfomin = 0.98
		self._lfomax = 1.02
		self._lfoseed = 1
		self._lfosmooth = 0
		self._lfodurrange = (0, 0)
		self._amplfotype = None
		self._amplfotype_israndom = False
		self._amplforate = 5.0
		self._amplfomin = 0.0
		self._amplfomax = 1.0
		self._amplfoseed = 1
		self._amplfosmooth = 0
		self._amplfodurrange = (0, 0)
		self._partials = []
		self._max_partial_index = 0
		self._num_unfiltered_partials = 0
		self._zero_phase = False;
		self._zero_first_partial_amp = False;
		self._zero_last_partial_amp = False;
		self._read_file_called = False
		self._pan_randgen = random.Random()
		self._pan_randgen.seed(1)
		self._synth_randgen = random.Random()
		self._synth_randgen.seed(2)
		self._delay_randgen = random.Random()
		self._delay_randgen.seed(3)
		self._quantize_randgen = random.Random()
		self._quantize_randgen.seed(4)
		self._gliss_randgen = random.Random()
		self._gliss_randgen.seed(5)
		self._lfo_randgen = random.Random()
		self._lfo_randgen.seed(6)
		self._amplfo_randgen = random.Random()
		self._amplfo_randgen.seed(7)
		rtcmix.load("WAVETABLE")
		rtcmix.load("BWESINE")

	def _is_number(self, thing):
		if type(thing) is types.IntType or type(thing) is types.FloatType:
			return True
		return False

	def _is_list(self, thing):
		if type(thing) is types.TupleType or type(thing) is types.ListType:
			return True
		return False

	def _is_rtcmix_handle(self, thing):
		if type(thing) is type(self._wavetable):
			return True
		return False

	""" Return a version string. """
	def version(self):
		return self._version

	""" If true, print out a bunch of stuff while parsing. """
	def set_verbose(self, verbose):
		self._verbose = verbose

	""" Exclude partials whose index numbers are not within this range.
	    Use max = 0 to place no upper limit on index.
	    Partial index has nothing to do with amp or frequency ordering.
	    Partials with greater indices tend to have later start times.
	    NB: Must call before read_file().
	"""
	def select_partial_range(self, min, max=0):
		if self._read_file_called:
			print "Must call select_partial_range before calling read_file."
			sys.exit()
		if max == 0 and self._num_unfiltered_partials != 0:
			max = self._num_unfiltered_partials
		self._partial_range = (min, max)

	""" Exclude partials whose start and end times are not within this range.
	    Use max = 0 to place no upper limit on time. If clipend is True, accept
		 a partial even if its end time is later than max, but clip the partial
		 so that it ends just before max.
	    NB: Must call before read_file().
	"""
	def select_time_range(self, min, max=0, clipend=False):
		if self._read_file_called:
			print "Must call select_time_range before calling read_file."
			sys.exit()
		self._time_range = (min, max)
		self._clip_end = clipend

	""" Exclude partials whose total durations are not all within this range.
	    Use max = 0 to place no upper limit on duration.
	    NB: Must call before read_file().
	"""
	def select_duration_range(self, min, max=0):
		if self._read_file_called:
			print "Must call select_duration_range before calling read_file."
			sys.exit()
		self._duration_range = (min, max)

	""" Exclude partials whose breakpoint amplitudes are not *all* within this
	    range. Use max = 0 to place no upper limit on amp. Keep in mind that
	    partial amplitudes are quite low (e.g., 0.00001 to 0.1). The final
	    breakpoint of a partial is often (always?) zero. We ignore this for the
	    purposes of determining exclusion.
	    NB: Must call before read_file().
	"""
	def select_amp_range(self, min, max=0):
		if self._read_file_called:
			print "Must call select_amp_range before calling read_file."
			sys.exit()
		self._amp_range = (min, max)

	""" Exclude partials whose breakpoint frequencies are not all within this
	    range. Use max = 0 to specify the Nyquist frequency as the upper limit.
	    NB: Must call before read_file().
	"""
	def select_freq_range(self, min, max=0):
		if self._read_file_called:
			print "Must call select_freq_range before calling read_file."
			sys.exit()
		if max == 0:
			max = self._nyquist
		self._freq_range = (min, max)

	""" Exclude partials whose frequency standard deviation is not within this
	    range. The std is a measure of how the individual breakpoint frequencies
	    deviate from the mean of all the partial frequencies. The larger the std,
	    the more pronounced is the glissando of the partial shape.
	    Use max = 0 to place no upper limit on std.
	    NB: Must call before read_file().
	"""
	def select_freq_std_range(self, min, max=0):
		if self._read_file_called:
			print "Must call select_freq_std_range before calling read_file."
			sys.exit()
		if max == 0:
			max = 1000
		self._freq_std_range = (min, max)

	""" Exclude partials that are not close enough to harmonic partials of the
	    given fundamental frequency (in Hz). Tolerance is a percentage in [0,100]
	    of the fundamental frequency; it determines how far away from a harmonic
	    partial a partial can be and still be retained. The determination is
	    made based on the mean partial frequency.
	    NB: Must call before read_file().
	"""
	def select_harmonic_partials(self, fundamental, tolerance=0):
		if self._read_file_called:
			print "Must call select_harmonic_partials before calling read_file."
			sys.exit()
		if fundamental < 20:
			print "WARNING: select_harmonic_partials: fundamental < 20 Hz"
			return
		harms = []
		partialnum = 1
		while True:
			thisfreq = fundamental * partialnum
			if thisfreq > self._nyquist:
				break
			harms.append(thisfreq)
			partialnum += 1
		self._freq_harmonics = harms
		self._freq_harm_tolerance = tolerance

	""" Exclude partials that are not close enough to any of the given chord
	    pitches.
	       chord        chord used for matching [pitches in oct.pc]
	       transp       transpose chord by this interval before matching
	                    [semitones]
	       sensitivity  retain partial if interval between mean partial
	                    frequency and nearest chord pitch is no greater
	                    than half this interval in either direction [semitones]
	    Note that this does not guarantee that you will hear chord notes when
	    excluding partials, because the retained partials may have high std
	    (the partial deviates a lot from its mean freq).
	    NB: Must call before read_file().
	"""
	def select_chord_partials(self, chord, transp=0, sensitivity=6):
		if self._read_file_called:
			print "Must call select_chord_partials before calling read_file."
			sys.exit()
		# cache chord data converted to linear octaves
		self._freq_chord_pitches = []
		for pitch in chord:
			self._freq_chord_pitches.append(rtcmix.octpch(pitch))
		self._freq_chord_transp = rtcmix.octpch(transp * 0.01)
		self._freq_chord_sensitivity = rtcmix.octpch(sensitivity * 0.01) * 0.5

	""" Exclude partials whose breakpoint bandwidths (if present) are not all
	    within this range. Note that a <max> of zero really means zero, not a
	    shorthand for the maximum permissable value, thus differing from the
	    behavior of parallel functions for duration, amplitude, and frequency.
	    This is to let you select only partials that have zero bandwidth.
	    NB: Must call before read_file().
	"""
	def select_bandwidth_range(self, min=0, max=1):
		if self._read_file_called:
			print "Must call select_bandwidth_range before calling read_file."
			sys.exit()
		self._bw_range = (min, max)

	def _do_read_file(self, filename):
		f = open(filename)
		lines = f.readlines()
		f.close()
		partials = []
		self._num_unfiltered_partials = totalpartials = 0
		pindex = numpoints = 0
		starttime = endtime = 0.0
		earliest_partial_start = sys.float_info.max
		numskipped = 0
		lineno = 1
		for line in lines:
			if lineno == 1:
				if line.startswith("par-text-partials-format"):
					self._extended = False
					print "Reading par-text-partials-format."
				elif line.startswith("par-text-partials-extended-format"):
					self._extended = True
					print "Reading par-text-partials-extended-format."
				else:
					print "Invalid file format type (line 1)."
			elif lineno == 2:
				# documentation of data items for each point
				pass
			elif lineno == 3:
				if not line.startswith("partials-count"):
					print "Invalid file format: no partials-count (line 3)."
				elems = line.split()
				totalpartials = elems[1]
			elif lineno == 4:
				pass
			else:
				elems = line.split()
				if lineno % 2:		# odd line number: header
					pindex = int(elems[0])		# numbers from file are strings
					self._max_partial_index = pindex
					numpoints = int(elems[1])
					starttime = float(elems[2])
					if starttime < 0.0:
						print "Negative starttime ({:.6f}) for index {}. Correcting.\n".format(starttime, pindex)
						starttime = 0.0
						# FIXME: but doesn't fix any breakpoints
					endtime = float(elems[3])
					#print "read: [{}] points={}, start={:.5}, end={:.5}".format(pindex, numpoints, starttime, endtime)
					if self._extended:
						phase = float(elems[4])
						label = int(elems[5])
					else:
						phase = 0
						label = 0
				else:					# even line number: partial data
					if self._partial_range[1] == 0:
						last_partial = totalpartials
					else:
						last_partial = self._partial_range[1]
					if pindex >= self._partial_range[0] and pindex <= last_partial:
						partial = []
						if self._extended:
							skip = 4
						else:
							skip = 3
						end = numpoints * skip
						for start in range(0, end, skip):
							time = float(elems[start])
							if time < 0.0:
								time = 0.0
							freq = float(elems[start + 1])
							amp = float(elems[start + 2])
							bp = []
							bp.append(time)
							bp.append(freq)
							bp.append(amp)
							if self._extended:
								bw = float(elems[start + 3])
							else:
								bw = 0
							bp.append(bw)
							partial.append(bp)
						if time != endtime:
							print "WARNING: Last breakpoint time is not the same as partial end time [partial {}].".format(pindex)
						partial = self._filter_partial(partial, starttime, endtime)
						if (partial == None):
							numskipped += 1
						else:
							endtime = partial[-1][0]	# may've been changed by _filter_partial
							mute = False
							p = [pindex, numpoints, starttime, endtime, phase, label, mute, partial]
							partials.append(p)
							self._num_unfiltered_partials += 1
							if starttime < earliest_partial_start:
								earliest_partial_start = starttime
			lineno += 1
		self._selected_partials_start = earliest_partial_start
		print "can play", self._num_unfiltered_partials, "out of", \
				totalpartials, "partials; skipping", numskipped

		# compute total input duration of job from start of first partial to
		# end of last partial
		# FIXME: is the last partial in list guaranteed to be temporally last?
		# and isn't this given in the header?
		if self._num_unfiltered_partials > 0:
			firstpartial = partials[0]
			numpartials = len(partials)
			lastpartial = partials[numpartials - 1]
			self._totdur = lastpartial[3] - firstpartial[2]
			print "total input duration:", self._totdur, "seconds"
		return partials

	""" Read a SPEAR par-text-partials-format file, store partials, and return
	    the number of unfiltered partials. Exit if trouble reading the file.
	"""
	def read_file(self, filename):
		self._read_file_called = True
		self._partials = self._do_read_file(filename)
		if len(self._partials) == 0:
			sys.exit()
		return self._num_unfiltered_partials

	""" Write a SPEAR par-text-partials-format file with the <partials> as the
	    data, or if that is None, the internal partials used by other methods,
	    including any transformations applied or specified earlier. Also supports
	    our par-text-partials-extended-format. Returns 0 if success, -1 if
	    failure.
	    NOTE1: Ignores a time-varying timescale.
	    NOTE2: Ignores prior calls to set_wavetable, set_ampenv2, 
	           set_glisstable, set_lfo, and set_amplfo.
	"""
	def write_file(self, filename, extended=False, partials=None):
		if partials == None:
			partials = self._partials
		lines = []
		header = partials[0]
		if extended:
			line = "par-text-partials-extended-format\n"
			lines.append(line)
			line = "point-type time frequency amplitude bandwidth\n"
			lines.append(line)
		else:
			line = "par-text-partials-format\n"
			lines.append(line)
			line = "point-type time frequency amplitude\n"
			lines.append(line)
		line = "partials-count {}\n".format(self._num_unfiltered_partials)
		lines.append(line)
		line = "partials-data\n"
		lines.append(line)
		index = 0
		for partial in partials:
			numpoints = partial[1]
			starttime = origstart = partial[2]
			endtime = partial[3]
			if self._timescale_start_points:
				starttime *= self._timescale
			endtime = starttime + ((endtime - origstart) * self._timescale)
			if extended:
				if self._zero_phase:
					phase = 0
				else:
					phase = partial[4]
				label = partial[5]
				line = "{} {} {:.6f} {:.6f} {:.6f} {}\n".format(index, numpoints, starttime, endtime, phase, label)
			else:
				line = "{} {} {:.6f} {:.6f}\n".format(index, numpoints, starttime, endtime)
			lines.append(line)
			breakpoints = partial[7]
			count = 1
			for bp in breakpoints:
				start = starttime + ((bp[0] - origstart) * self._timescale)
				freq = bp[1]
				if (count == 1 and self._zero_first_partial_amp) \
						or (count == numpoints and self._zero_last_partial_amp):
					amp = 0.0
				else:
					amp = bp[2] * self._ampscale
				if extended:
					bw = bp[3]
					line = "{:.6f} {:.6f} {:.6f} {:.6f} ".format(start, freq, amp, bw)
				else:
					line = "{:.6f} {:.6f} {:.6f} ".format(start, freq, amp)
				lines.append(line)
				count += 1
			lines[-1].rstrip(' ')	# chomp terminating space
			lines.append('\n')
			index += 1
		if index != self._num_unfiltered_partials:
			print "WARNING: written partial index does not match number of unfiltered partials"
		f = open(filename, "w")
		f.writelines(lines)
		f.close()
		return 0

	""" Helper for print_histogram. Return index into octave band, with
	    0 for anything below 22 Hz, the bottom edge of the first band;
	    1 for anything within the 31.5 Hz band (22-44 Hz).
	    Band edges derived from
	    https://engineeringtoolbox.com/octave-bands-frequency-limits-d_1602.html
	    Bands confirmed by https://en.wikipedia.org/wiki/Octave_band.
	    Exception: our first band extends from 0 Hz, not 11 Hz.
	"""
	def _octave_band_index(self, freq):
		#         16 31.5 63  125 250  500  1000  2000  4000  8000  16000
		bands = (0, 22, 44, 88, 177, 355, 710, 1420, 2840, 5680, 11360, 99999)
		for i in range(0, len(bands)):
			if bands[i] > freq:
				return i - 1
		return -1	# can't happen until srates go above 192 kHz

	""" Print to the screen a histogram that shows the number of partials for
	    each octave band having amplitudes, durations, bandwidths, and frequency
	    standard deviations within the given ranges. Here is the definition of
	    the frequency bands:
	          16 Hz:     0 - 22
	          31 Hz:    22 - 44
	          63 Hz:    44 - 88
	         125 Hz:    88 - 177
	         250 Hz:   177 - 355
	         500 Hz:   355 - 710
	        1000 Hz:   710 - 1420
	        2000 Hz:  1420 - 2840
	        4000 Hz:  2840 - 5680
	        8000 Hz:  5680 - 11360
	       16000 Hz: 11360 - 99999
	    The histogram includes only partials that have not been filtered out.
	    Partial membership in a frequency band is determined by the mean
	    frequency of the partial, so it might drift outside of this band.
	"""
	def print_histogram(self, minamp, maxamp, mindur, maxdur, minbw, maxbw, minstd, maxstd, partials=None):
		if partials == None:
			partials = self._partials
		bandnames = ("16", "31", "63", "125", "250", "500", "1000", "2000", \
						"4000", "8000", "16000")
		ampcounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		durcounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		bwcounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		stdcounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		tminamp = tmindur = tminbw = tminstd = sys.float_info.max
		tmaxamp = tmaxdur = tmaxbw = tmaxstd = -1.0
		for p in partials:
			starttime = p[2]
			endtime = p[3]
			breakpoints = p[7]
			meanfreq = self.get_partial_freq_mean(breakpoints)
			band = self._octave_band_index(meanfreq)
			# amp
			(pminamp, pmaxamp) = self.get_partial_amp_range(breakpoints)
			if pminamp >= minamp and pmaxamp <= maxamp:
				ampcounts[band] += 1
			if tminamp > pminamp:
				tminamp = pminamp
			if tmaxamp < pmaxamp:
				tmaxamp = pmaxamp
			# dur
			dur = endtime - starttime
			if dur >= mindur and dur <= maxdur:
				durcounts[band] += 1
			if tmindur > dur:
				tmindur = dur
			if tmaxdur < dur:
				tmaxdur = dur
			# bw
			(pminbw, pmaxbw) = self.get_partial_bw_range(breakpoints)
			if pminbw >= minbw and pmaxbw <= maxbw:
				bwcounts[band] += 1
			if tminbw > pminbw:
				tminbw = pminbw
			if tmaxbw < pmaxbw:
				tmaxbw = pmaxbw
			# std
			std = self.get_partial_freq_std(breakpoints, meanfreq)
			if std >= minstd and std <= maxstd:
				stdcounts[band] += 1
			if tminstd > std:
				tminstd = std
			if tmaxstd < std:
				tmaxstd = std
		# print results
		print "\nHistogram: num. partials with values in requested ranges"
		print   "---- requested ----------\tactual -------------"
		print "amp: [{:f}, {:f}]\t[{:f}, {:f}]".format(minamp, maxamp, tminamp, tmaxamp)
		print "dur: [{:f}, {:f}]\t[{:f}, {:f}]".format(mindur, maxdur, tmindur, tmaxdur)
		print " bw: [{:f}, {:f}]\t[{:f}, {:f}]".format(minbw, maxbw, tminbw, tmaxbw)
		print "std: [{:f}, {:f}]\t[{:f}, {:f}]".format(minstd, maxstd, tminstd, tmaxstd)
		print "------------------------------------------------"
		print "    band       amp       dur        bw       std"
		print "------------------------------------------------"
		for i in range(0, len(ampcounts)):
			print "{:>5} Hz:   {:6d}    {:6d}    {:6d}    {:6d}".format(bandnames[i], ampcounts[i], durcounts[i], bwcounts[i], stdcounts[i])
		print ""

	""" Print to the screen the parameter ranges exhibited by partials within
	    the specified frequency range. The time range and the collection of
	    partials under consideration are affected by whatever select methods
	    already have been called prior to printing the frequency band.
	    Note that a partial must be entirely contained within the current
	    time range to be included here.
	"""
	def print_freq_range_params(self, minfreq, maxfreq=0, partials=None):
		if partials == None:
			partials = self._partials
		if maxfreq == 0:
			maxfreq = self._nyquist
		req_minfreq = minfreq
		req_maxfreq = maxfreq
		req_minstart = self._time_range[0]
		req_maxend = self._time_range[1]
		minfreq = sys.float_info.max
		maxfreq = 0
		minstart = sys.float_info.max
		maxend = 0
		minamp = sys.float_info.max
		maxamp = 0
		mindur = sys.float_info.max
		maxdur = 0
		minbw = sys.float_info.max
		maxbw = 0
		minstd = sys.float_info.max
		maxstd = 0
		count = 0
		for p in partials:
			breakpoints = p[7]
			# freq
			(pminfreq, pmaxfreq) = self.get_partial_freq_range(breakpoints)
			if minfreq > pminfreq:
				minfreq = pminfreq
			if maxfreq < pmaxfreq:
				maxfreq = pmaxfreq
			if pminfreq < req_minfreq or pmaxfreq > req_maxfreq:
				continue
			# time
			starttime = p[2]
			endtime = p[3]
			if starttime < minstart:
				minstart = starttime
			if endtime > maxend:
				maxend = endtime
			# amp
			(pminamp, pmaxamp) = self.get_partial_amp_range(breakpoints)
			if minamp > pminamp:
				minamp = pminamp
			if maxamp < pmaxamp:
				maxamp = pmaxamp
			# dur
			dur = endtime - starttime
			if mindur > dur:
				mindur = dur
			if maxdur < dur:
				maxdur = dur
			# bw
			(pminbw, pmaxbw) = self.get_partial_bw_range(breakpoints)
			if minbw > pminbw:
				minbw = pminbw
			if maxbw < pmaxbw:
				maxbw = pmaxbw
			# std
			meanfreq = self.get_partial_freq_mean(breakpoints)
			std = self.get_partial_freq_std(breakpoints, meanfreq)
			if minstd > std:
				minstd = std
			if maxstd < std:
				maxstd = std
			count += 1
		if count == 0:
			minstart = req_minstart
			maxend = req_maxend
		# print results
		print "\nPrint params for time range: {:.3f} - {:.3f}".format(req_minstart, req_maxend)
		print   "             and freq range: {:.3f} - {:.3f}".format(req_minfreq, req_maxfreq)
		print   "{:d} partials selected for printing, in these ranges...".format(count)
		print "---------------------------------------------------"
		print "time: [{:f}, {:f}]".format(minstart, maxend)
		print "freq: [{:f}, {:f}]".format(minfreq, maxfreq)
		if count > 0:
			print " amp: [{:f}, {:f}]".format(minamp, maxamp)
			print " dur: [{:f}, {:f}]".format(mindur, maxdur)
			print "  bw: [{:f}, {:f}]".format(minbw, maxbw)
			print " std: [{:f}, {:f}]".format(minstd, maxstd)
		print "---------------------------------------------------"
		print ""

	""" Use this is you have recompiled RTcmix with a greater number for
	    MAXDISPARGS in rtcmix/src/include/maxdispargs.h. This lets you work
	    with partials containing a larger number of breakpoints. The default
	    is 1024, which gives you 510 breakpoints. Use a power of two for this.
	"""
	def set_max_disp_args(self, max=1024):
		# time,val pairs, and leaving room for table type, "nonorm", size, etc.
		self._max_num_breakpoints = (max - 4) / 2

	""" Set the RTcmix output bus for Resynth. Do this only for mono and stereo.
	    Multichannel output works only by choosing point-source output channels
	    randomly per partial. For this, use set_multichan_bustype.
	"""
	def set_outbus(self, bus="out 0-1"):
		# do both, because we use WAVETABLE if a partial has zero noise bandwidth
		rtcmix.bus_config("WAVETABLE", bus)
		rtcmix.bus_config("BWESINE", bus)

	""" If the number of output channels is greater than 2, use this bus type
	    for the RTcmix instruments playing the partials. Each partial will come
	    out of a single randomly-selected bus of the type given, which can be
	    either "out" (default) or "aux". Bus numbers are 0 through numchans-1,
	    where numchans is the second argument to the Resynth constructor.
	"""
	def set_multichan_bustype(self, bustype="out"):
		if bustype is "out" or bustype is "aux":
			self._multichan_bustype = bustype
		else:
			print "\nERROR: set_multichan_bustype: invalid bus type '{:s}'.".format(bustype)
			sys.exit()

	""" Set overall volume level in dB (0 = no change). """
	def set_gain(self, gain=0):
		self._ampscale = rtcmix.ampdb(gain)

	""" Pin amplitudes to the given range, whose min and max cannot be negative.
	    Amplitudes already within this range are not altered.
	"""
	def clamp_amps(self, min, max, partials=None):
		if partials == None:
			partials = self._partials
		if min < 0.0:
			print "clamp_amps: min must be => 0.0."
			return partials
		if max < 0.0:
			print "clamp_amps: max must be => 0.0."
			return partials
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				if bp[2] < min:
					bp[2] = min
				elif bp[2] > max:
					bp[2] = max
		return partials

	""" Invert amplitudes with respect to the max amp value for the entire
	    set of selected partials.
	"""
	def invert_amps(self, partials=None):
		if partials == None:
			partials = self._partials
		max = -sys.float_info.max
		for p in partials:
			breakpoints = p[7]
			(thismin, thismax) = self.get_partial_amp_range(breakpoints)
			if thismax > max:
				max = thismax
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				bp[2] = max - bp[2]
		return partials

	""" Scale noise bandwidth of all partials. This has an effect only if
	    reading BWE partials in our extended partext format. Note that the
	    RTcmix BWESINE instrument clamps the bandwidth to the range [0,1].
	    The scale argument can be either a constant or an RTcmix PFIeld
	    table handle, where the table spans the entire original duration
	    of the selected set of partials.
	"""
	def scale_bandwidths(self, scale, partials=None):
		if partials == None:
			partials = self._partials
		tablen = 0
		if type(scale) is not int and type(scale) is not float:
			# we assume it's an RTcmix-specific value
			isTable = True
			tablen = rtcmix.tablelen(scale)
		else:
			isTable = False
			if scale < 0.0:
				print "scale_bandwidths: scale must not be less than 0.0."
				return partials
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				if isTable:
					start = bp[0] - self._selected_partials_start
					index = (start / self._totdur) * tablen
					bp[3] *= rtcmix.samptable(scale, index)
				else:
					bp[3] *= scale
		return partials

	""" Pin noise bandwidths to the given range, which must be between
	    0.0 and 1.0. Bandwidths already within this range are not altered.
	"""
	def clamp_bandwidths(self, min, max, partials=None):
		if partials == None:
			partials = self._partials
		if min < 0.0 or min > 1.0:
			print "clamp_bandwidths: min must be between 0.0 and 1.0."
			return partials
		if max < 0.0 or max > 1.0:
			print "clamp_bandwidths: max must be between 0.0 and 1.0."
			return partials
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				if bp[3] < min:
					bp[3] = min
				elif bp[3] > max:
					bp[3] = max
		return partials

	""" Shift the time points for all breakpoints by <shift>.
	    Note that partials with a negative start time will not be played.
	"""
	def shift_times(self, shift, partials=None):
		if partials == None:
			partials = self._partials
		for p in partials:
			p[2] += shift	# starttime
			p[3] += shift	# endtime
			breakpoints = p[7]
			for bp in breakpoints:
				bp[0] += shift
		self._selected_partials_start += shift
		return partials

	""" Delay the time points for all breakpoints in the given frequency range.
	       delay        an RTcmix table handle, with table values as delay times
	                    (in seconds) and table indices as frequency from 0 Hz to
	                    Nyquist. The frequencies are the amp-weighted mean
	                    frequencies of partials. Use the "nonorm" flag when
	                    creating this table.
	       deldev       maximum random deviation (in seconds) applied to a
	                    breakpoint time after delaying it (seeded using
	                    set_delay_seed method)
	       fmin, fmax   range of frequencies to affect (fmax=0 means Nyquist);
	                    a partial is delayed if its amp-weighted mean frequency
	                    lies within this range (inclusive).
	       partials     list of partials to delay [default is all]
	"""
	def delay_times(self, delay, deldev=0, fmin=0, fmax=0, partials=None):
		if partials == None:
			partials = self._partials
		if fmax == 0:
			fmax = self._nyquist
		if not self._is_rtcmix_handle(delay):
			print "delay_times: delay must be an RTcmix table handle."
			sys.exit()
		index_factor = rtcmix.tablelen(delay) / self._nyquist
		for p in partials:
			breakpoints = p[7]
			mean = self.get_partial_freq_mean(breakpoints)
			if mean >= fmin and mean <= fmax:
				shift = rtcmix.samptable(delay, mean * index_factor)
				if deldev != 0:
					dev = self._delay_randgen.uniform(-deldev, deldev)
					shift += dev
				p[2] += shift	# starttime
				p[3] += shift	# endtime
				for bp in breakpoints:
					bp[0] += shift
		return partials

	""" Quantize the start time for every partial in the given frequency range,
	    and shift all of its breakpoint times accordingly.
	       quantum      partial start times will be integer multiples of this
	       quantdev     maximum random deviation (in seconds) applied to a
	                    partial start time after quantizing it (seeded using
	                    set_quantize_seed method). This affects equally all
	                    partials that quantize to the same attack point.
	       smear        additional random deviation (in seconds) applied to 
	                    individual partials that find themselves starting at the
	                    same attack point (seeded using set_quantize_seed method)
	       template     optional list or tuple of numbers that determine whether
	                    the partials at a specific quantized startpoint will
	                    play. Think of the sequence of numbers as aligned with
	                    the sequence of quantized attack points (i.e., integer
	                    multiples of the quantum value). If a given template
	                    value is zero, the partials beginning at that start point
	                    will not play; otherwise they will. The template sequence
	                    repeats as long as necessary to span the series of
	                    quantized start points. The template can be any length.
	       fmin, fmax   range of frequencies to affect (fmax=0 means Nyquist);
	                    a partial is quantized if its amp-weighted mean frequency
	                    lies within this range (inclusive).
	       partials     list of partials to quantize [default is all]
	"""
	def quantize_times(self, quantum, quantdev=0, smear=0, template=0, fmin=0, fmax=0, partials=None):
		if partials == None:
			partials = self._partials
		if fmax == 0:
			fmax = self._nyquist
		# Cache each unique quantized start point, and associate with it one
		# randomized deviation value. Apply this deviation to all partials
		# that are quantized to that start point. Then apply smear deviation
		# to all partials, independently of their start points.
		starttime_cache = {}
		precision = 6
		for p in partials:
			breakpoints = p[7]
			mean = self.get_partial_freq_mean(breakpoints)
			if mean >= fmin and mean <= fmax:
				origstart = float(p[2])
				quotient = abs(origstart / quantum)
				floor = int(quotient)
				remainder = quotient - floor
				if remainder >= 0.5:
					floor += 1.0
				newstart = floor * quantum
				if origstart < 0.0:
					newstart *= -1.0
				shift = newstart - origstart
				if quantdev != 0:
					key = round(newstart, precision) # hash consistency for floats
					if key in starttime_cache:
						dev = starttime_cache[key]
					else:
						dev = self._quantize_randgen.uniform(-quantdev, quantdev)
						starttime_cache[key] = dev
					shift += dev
				if smear != 0:
					dev = self._quantize_randgen.uniform(-smear, smear)
					shift += dev
				if template != 0:
					if not self._is_list(template):
						print "quantize_times: template must be 0 or be a list or tuple of numbers."
						sys.exit()
					index = int(newstart / quantum) % len(template)
					if template[index] == 0:
						p[6] = True   # mute
					else:
						p[6] = False  # reset
				p[2] += shift	# starttime
				p[3] += shift	# endtime
				for bp in breakpoints:
					bp[0] += shift
		return partials

	""" If true, set all partial starting phases to zero (bool). """
	def set_zero_phase(self, doit):
		self._zero_phase = doit

	""" If true, zero out the amplitude of the first breakpoint of every partial;
	    otherwise, leave the amp as it was in the SPEAR file.
	"""
	def set_zero_first_partial_amp(self, doit):
		self._zero_first_partial_amp = doit;

	""" If true, zero out the amplitude of the final breakpoint of every partial;
	    otherwise, leave the amp as it was in the SPEAR file. Often this is zero,
	    but not always.
	"""
	def set_zero_last_partial_amp(self, doit):
		self._zero_last_partial_amp = doit;

	""" Stretch partial durations by this factor, which can be a constant or
	    an RTcmix table handle.
	    Behavior depends on set_timescale_start_points().
	"""
	def set_timescale(self, timescale):
		if self._is_number(timescale):
			self._timescale = timescale
		else:
			print "Resynth does not yet support time-varying timescale."
			sys.exit()

	""" If false, don't change partial start points when time-scaling;
	    else, scale start points (bool). Default is true.
	"""
	def set_timescale_start_points(self, doit):
		self._timescale_start_points = doit

	""" Set an RTcmix wavetable to use in place of a sine wave. """
	def set_wavetable(self, table_handle):
		self._wavetable = table_handle

	""" Set a curve to scale the amp envelope of each partial. This can be 
	    either an RTcmix table handle, whose duration will be scaled to fit
	    the duration of each partial, or a tuple of four arguments that 
	    help generate an RTcmix "curve" table for each partial.

	       attacktime      time to reach full amplitude breakpoint
	       attackcurve     shape of curve approaching this breakpoint
	       holdtime        time to hold at max amplitude before release
	       releasecurve    curvature down to zero amplitude

	    The advantage of this way is that the times above will be constant
	    across partials, regardless of their durations. The exception is
	    that these times will be shortened when the partial is too short to
	    accommodate them. Specifically, the attack time is preserved as far
	    as possible, sacrificing hold time, but without producing a release
	    click. The curve amounts are in the same format as those given to
	    the "curve" table in RTcmix.
	"""
	def set_ampenv2(self, table):
		islist = self._is_list(table)
		if islist and len(table) != 4:
			print "set_ampenv2 table spec takes 4 arguments."
			sys.exit()
		if islist or self._is_rtcmix_handle(table):
			self._ampenv2 = table
		else:
			print "set_ampenv2: table must be a valid list or an RTcmix table."
			sys.exit()

	def set_eqtable(self, table_handle):
		print "set_eqtable is deprecated; use eq_partials instead."

	""" Set an RTcmix table (or tables -- see below) to use for time-varying
	    scaling of the frequencies of a partial. The table values are frequency
	    multipliers, with the duration of the table scaled to fit the partial
	    duration. Use the "nonorm" flag to avoid accidental multipliers of 0.
	    Applies the table to all partials whose durations, post timescale, are
	    between mindur and maxdur (where maxdur=0 specifies no upper boundary).
	    If table_handle is a list or tuple of table handles, then these are
	    randomly selected during playback. Use set_gliss_seed to change the
	    sequence of tables selected.
	"""
	def set_glisstable(self, table_handle, mindur=0, maxdur=0):
		self._glisstable = table_handle
		self._glissdurrange = (mindur, maxdur)

	""" Specify an LFO to create for the frequency of each partial, using values
	    that will be passed to RTcmix makeLFO or makerandom, depending on the
	    supplied type string (which is the waveform string or handle for makeLFO
	    and random distribution type string for makerandom). The only type not
	    supported is makerandom "prob". <rate> is either a constant (Hz) or an
	    RTcmix table handle. <min> and <max> are multipliers that will affect
	    the time-varying frequency of the partial. <rate>, <min>, and <max> can
	    be tuples or lists of constants, as well as single constants. Use
	    set_lfo_seed to change the sequence of values pulled from these
	    tuples/lists. <rate>, <min>, and <max> can also be RTcmix table handles.
	    <seed> and <smooth> are used only for the random types: <seed> is passed
	    to rtcmix.makerandom, while <smooth> is a percentage for the smoothing
	    filter that processes the resulting random number stream. The LFO
	    applies to all partials whose durations, post timescale, are between
	    mindur and maxdur (where maxdur=0 specifies no upper boundary). The LFO
	    is a relatively expensive capability.
	"""
	def set_lfo(self, type, rate, min, max, seed=1, smooth=0, mindur=0, maxdur=0):
		self._lfotype = type
		if self._is_rtcmix_handle(type):		# waveform in a table handle
			self._lfotype_israndom = False
		elif type is "even" or type is "linear" or type is "low" \
				or type is "high" or type is "triangle" or type is "gaussian" \
				or type is "cauchy":
			self._lfotype_israndom = True
		else:
			self._lfotype_israndom = False
		self._lforate = rate
		self._lfomin = min
		self._lfomax = max
		self._lfoseed = seed
		self._lfosmooth = smooth
		self._lfodurrange = (mindur, maxdur)

	""" Specify an LFO to create for the amplitude of each partial, using values
	    that will be passed to RTcmix makeLFO or makerandom, depending on the
	    supplied type string (which is the waveform string or handle for makeLFO
	    and random distribution type string for makerandom). The only type not
	    supported is makerandom "prob". <rate> is either a constant (Hz) or an
	    RTcmix table handle. <min> and <max> are multipliers that will affect
	    the time-varying amplitude of the partial. <rate>, <min>, and <max> can
	    be tuples or lists of constants, as well as single constants. Use
	    set_amplfo_seed to change the sequence of values pulled from these
	    tuples/lists. <rate>, <min>, and <max> can also be RTcmix table handles.
	    <seed> and <smooth> are used only for the random types: <seed> is passed
	    to rtcmix.makerandom, while <smooth> is a percentage for the smoothing
	    filter that processes the resulting random number stream. The LFO
	    applies to all partials whose durations, post timescale, are between
	    mindur and maxdur (where maxdur=0 specifies no upper boundary). The LFO
	    is a relatively expensive capability.
	"""
	def set_amplfo(self, type, rate, min, max, seed=1, smooth=0, mindur=0, maxdur=0):
		self._amplfotype = type
		if self._is_rtcmix_handle(type):		# waveform in a table handle
			self._amplfotype_israndom = False
		elif type is "even" or type is "linear" or type is "low" \
				or type is "high" or type is "triangle" or type is "gaussian" \
				or type is "cauchy":
			self._amplfotype_israndom = True
		else:
			self._amplfotype_israndom = False
		self._amplforate = rate
		self._amplfomin = min
		self._amplfomax = max
		self._amplfoseed = seed
		self._amplfosmooth = smooth
		self._amplfodurrange = (mindur, maxdur)

	""" Seed the random number generator used for panning when using stereo
	    output. Without this call the seed will be set to 1.
	"""
	def set_pan_seed(self, seed):
		self._pan_randgen.seed(seed)

	""" Seed the random number generator used for warping frequencies in the
	    synthesize_new_partials method. Without this call the seed will be
	    set to 2.
	"""
	def set_synth_seed(self, seed):
		self._synth_randgen.seed(seed)

	""" Seed the random number generator used for delaying partials in the
	    delay_times method. Without this call the seed will be set to 3.
	"""
	def set_delay_seed(self, seed):
		self._delay_randgen.seed(seed)

	""" Seed the random number generator used for quantizing partials in the
	    quantize_times method. Without this call the seed will be set to 4.
	"""
	def set_quantize_seed(self, seed):
		self._quantize_randgen.seed(seed)

	""" Seed the random number generator used for selecting from multiple
	    gliss tables. Without this call the seed will be set to 5.
	"""
	def set_gliss_seed(self, seed):
		self._gliss_randgen.seed(seed)

	""" Seed the random number generator used for selecting from multiple
	    frequency LFO values for rate, min, and max. Without this call the
	    seed will be set to 6.
	"""
	def set_lfo_seed(self, seed):
		self._lfo_randgen.seed(seed)

	""" Seed the random number generator used for selecting from multiple
	    amplitude LFO values for rate, min, and max. Without this call the
	    seed will be set to 7.
	"""
	def set_amplfo_seed(self, seed):
		self._amplfo_randgen.seed(seed)

	""" Get the list of partials, so that you can mangle it
	    and then pass it to modify_all_partials.
	"""
	def get_partials(self):
		return self._partials

	""" Return the time boundaries of the set of partials. """
	def get_time_bounds(self, partials=None):
		if partials == None:
			partials = self._partials
		minstart = sys.float_info.max
		maxend = -sys.float_info.max
		for p in partials:
			starttime = p[2]
			endtime = p[3]
			if starttime < minstart:
				minstart = starttime
			if endtime > maxend:
				maxend = endtime
		return (minstart, maxend)

	""" Return the freq range for the given set of partial breakpoints. """
	def get_partial_freq_range(self, partial_breakpoints):
		minfreq = sys.float_info.max
		maxfreq = 0.0
		for bp in partial_breakpoints:
			freq = bp[1]
			if freq >= 0.0 and freq < minfreq:
				minfreq = freq
			if freq > maxfreq:
				maxfreq = freq
		return (minfreq, maxfreq)

	""" Return the average frequency of a partial, weighted by the
	    amplitudes associated with each partial breakpoint frequency.
	    If <weighted> is false, return unweighted mean.
	"""
	def get_partial_freq_mean(self, partial_breakpoints, weighted=True):
		if weighted:
			freqs = 0.0
			amps = 0.0
			for bp in partial_breakpoints:
				freqs += bp[1] * bp[2]
				amps += bp[2]
			if amps == 0.0:
				# assuming positive amps, all freqs have 0 weight; avoid div/0 below
				return freqs / len(partial_breakpoints)
			return freqs / amps
		else:
			total = 0.0
			for bp in partial_breakpoints:
				total += bp[1]
			return total / len(partial_breakpoints)

	""" Return the standard deviation of the partial frequencies. """
	def get_partial_freq_std(self, partial_breakpoints, mean=None):
		if mean == None:	# if not already calculated
			mean = self.get_partial_freq_mean(partial_breakpoints)
		diffs = 0.0
		for bp in partial_breakpoints:
			diff = bp[1] - mean
			diffs += pow(diff, 2)
		variance =  diffs / len(partial_breakpoints)
		return math.sqrt(variance)	# std

	""" Return the amp range for the given set of partial breakpoints. """
	def get_partial_amp_range(self, partial_breakpoints):
		minamp = sys.float_info.max
		maxamp = 0.0
		for bp in partial_breakpoints:
			amp = bp[2]
			if amp > 0.0 and amp < minamp:
				minamp = amp
			if amp > maxamp:
				maxamp = amp
		return (minamp, maxamp)

	""" Return the bw range for the given set of partial breakpoints. """
	def get_partial_bw_range(self, partial_breakpoints):
		minbw = sys.float_info.max
		maxbw = 0.0
		for bp in partial_breakpoints:
			bw = bp[3]
			if bw < minbw:
				minbw = bw
			if bw > maxbw:
				maxbw = bw
		return (minbw, maxbw)

	""" Retune a single partial. """
	def _retune_one_partial(self, partial_breakpoints):
		# Is there a pitch in the chord list that is within _retune_sensitivity
		# semitones of the average pitch of this partial? NB: all pitches as
		# linear octaves, so that we can perform arithmetic on them.
		srcpitch = rtcmix.octcps(self.get_partial_freq_mean(partial_breakpoints))
		numchordpitches = len(self._retune_chord)
		mininterval = 999999.0
		index = -1
		for i in range(numchordpitches):
			targetpitch = self._retune_chord[i] + self._retune_chord_transp
			interval = abs(srcpitch - targetpitch)
			if (interval < mininterval):
				mininterval = interval
				index = i
		if (index >= 0 and mininterval <= self._retune_sensitivity):
			# if so, retune the partial
			if self._is_rtcmix_handle(self._retune_strength): # assume table
				tablen = rtcmix.tablelen(self._retune_strength)
			targetpitch = self._retune_chord[index] + self._retune_chord_transp
			# move each freq in the partial breakpoint list closer to the target
			# pitch, according to _retune_strength
			for bp in partial_breakpoints:
				if self._is_number(self._retune_strength):
					strength = self._retune_strength
				else: # assume table
					start = bp[0] - self._selected_partials_start
					index = (start / self._totdur) * tablen
					strength = rtcmix.samptable(self._retune_strength, index)
				origpitch = rtcmix.octcps(bp[1])
				diff = targetpitch - origpitch
				newpitch = origpitch + (diff * strength)
				bp[1] = rtcmix.cpsoct(newpitch) # retuned freq
			self._numretunedpartials += 1
		return partial_breakpoints

	""" Retune partials to the given chord.
	       chord        chord used for retuning [pitches in oct.pc]
	       transp       transpose chord by this interval [semitones]
	       sensitivity  retune if interval between partial frequency and
	                    nearest chord pitch is no greater than half of this
	                    interval in either direction [semitones]
	       strength     the extent to which a retuned partial conforms to a
                       close target pitch in the retune chord [0-1 scalar, or
                       RTcmix table handle]
	       partials     list of partials to retune [default is all]
	"""
	def retune_partials(self, chord, transp=0, sensitivity=6, strength=0.9, partials=None):
		# cache chord data converted to linear octaves
		self._retune_chord = []
		for pitch in chord:
			self._retune_chord.append(rtcmix.octpch(pitch))
		self._retune_chord_transp = rtcmix.octpch(transp * 0.01)
		self._retune_sensitivity = rtcmix.octpch(sensitivity * 0.01) * 0.5
		if self._is_number(strength):
			if (strength < 0 or strength > 1):
				print "retune strength must be between 0 and 1"
				sys.exit()
		self._retune_strength = strength
		if partials == None:
			partials = self._partials
		for p in partials:
			breakpoints = p[7]
			breakpoints = self._retune_one_partial(breakpoints)
		print self._numretunedpartials, "retuned partials"
		return partials

	""" Scale, then offset, all partial frequencies.
	       scale        multiply partial frequencies by this factor *
	       offset       then add this offset to partial frequencies *
	       freqdev      maximum pitch shift to apply randomly, up or down, to
	                    each scaled and offset frequency value [oct.pc semitones]
	       fmin, fmax   range of frequencies to affect (fmax=0 means Nyquist)
	       partials     list of partials to scale and offset [default is all]

	       * <scale> and <offset> can be RTcmix tables, with table values as
	         scale or offset and table indices as frequency from 0 Hz to Nyquist.
	         This lets the frequency warping affect different frequency ranges
	         differently. Use the "nonorm" flag when creating these tables.

	       Every breakpoint frequency is evaluated individually, without concern
	       for the partial it's in (e.g., the partial mean freq). Resulting
	       frequencies are limited to 20Hz - Nyquist.
	"""
	def scale_and_offset_freqs(self, scale, offset=0, freqdev=0, fmin=0, fmax=0, partials=None):
		if partials == None:
			partials = self._partials
		if fmax == 0:
			fmax = self._nyquist
		scale_is_table = False
		if self._is_rtcmix_handle(scale):
			scale_is_table = True
			scale_index_factor = rtcmix.tablelen(scale) / self._nyquist
		offset_is_table = False
		if self._is_rtcmix_handle(offset):
			offset_is_table = True
			offset_index_factor = rtcmix.tablelen(offset) / self._nyquist
		freqdev = abs(rtcmix.octpch(freqdev))
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				if bp[1] >= fmin and bp[1] <= fmax:
					if scale_is_table:
						s = rtcmix.samptable(scale, bp[1] * scale_index_factor)
					else:
						s = scale
					if offset_is_table:
						o = rtcmix.samptable(offset, bp[1] * offset_index_factor)
					else:
						o = offset
					bp[1] *= s
					bp[1] += o
					if freqdev != 0.0:
						target = rtcmix.octcps(bp[1])
						v = self._synth_randgen.uniform(-freqdev, freqdev)
						bp[1] = rtcmix.cpsoct(target + v)
					if bp[1] < 20:
						bp[1] = 20
					elif bp[1] > self._nyquist:
						bp[1] = self._nyquist
		return partials

	""" EQ partial breakpoints within the frequency range [fmin, fmax] using the
	    supplied RTcmix table handle. The table values are in dB and are
	    deployed along a linear frequency axis, with the first value assumed to
	    affect 0 Hz and the last value assumed to affect Nyquist. The recommended
	    approach is to use a "line" or "curve" table with 0 and Nyquist as the
	    first and last "time" values of the table parameters. Use the "nonorm"
	    flag in order to specify dB values with 0 meaning no change.

	    NOTE: Every breakpoint frequency is evaluated individually, without
	          concern for the partial it's in (e.g., the partial mean freq).
	"""
	def eq_partials(self, eqtable, fmin=0, fmax=0, partials=None):
		if partials == None:
			partials = self._partials
		if fmax == 0:
			fmax = self._nyquist
		if not self._is_rtcmix_handle(eqtable):
			print "eq_partials: eqtable must be an RTcmix table handle."
			sys.exit()
		index_factor = rtcmix.tablelen(eqtable) / self._nyquist
		for p in partials:
			breakpoints = p[7]
			for bp in breakpoints:
				if bp[1] >= fmin and bp[1] <= fmax:
					gain = rtcmix.samptable(eqtable, bp[1] * index_factor)
					bp[2] *= rtcmix.ampdb(gain)
		return partials

	""" For every partial whose mean freq is in the range defined by srcminfreq
	    and srcmaxfreq, synthesize a new partial whose frequencies are freqscale
	    multiples of the original freqs. Also apply gain (in dB) to each new
	    breakpoint. <freqdev> (in oct.pc semitones) is the maximum pitch shift
	    to apply randomly, up or down, to the transposed frequency. <minfreq>
	    is the lowest allowable frequency; high frequencies are limited to
	    Nyquist.
	"""
	def synthesize_new_partials(self, srcminfreq, srcmaxfreq, freqscale, gain=0, freqdev=0, minfreq=20, partials=None):
		if partials == None:
			partials = self._partials
		freqdev = abs(rtcmix.octpch(freqdev))
		ampscale = rtcmix.ampdb(gain)
		tmppartials = []
		for p in partials:
			breakpoints = p[7]
			mean = self.get_partial_freq_mean(breakpoints)
			if mean >= srcminfreq and mean <= srcmaxfreq:
				newpartial = copy.deepcopy(p)
				breakpoints = newpartial[7]
				for bp in breakpoints:
					bp[1] = bp[1] * freqscale
					if freqdev != 0.0:
						target = rtcmix.octcps(bp[1])
						v = self._synth_randgen.uniform(-freqdev, freqdev)
						bp[1] = rtcmix.cpsoct(target + v)
					if bp[1] < minfreq:
						bp[1] = minfreq
					elif bp[1] > self._nyquist:
						bp[1] = self._nyquist
					bp[2] = bp[2] * ampscale
				self._max_partial_index += 1
				newpartial[0] = self._max_partial_index
				tmppartials.append(newpartial)
		for partial in tmppartials:
			self._partials.append(partial)
			self._num_unfiltered_partials += 1
		return partials

	""" Scale the time points of all partial breakpoints, relative to the
	    start time for the partial, which is left unaltered. Partials will
	    still be subject to any time-scaling configured by set_timescale.
	"""
	def scale_time_points(self, scale, partials=None):
		if partials == None:
			partials = self._partials
		if scale <= 0.0:
			print "scale_time_points: scale must be greater than zero."
			return partials
		#print "----------------------------- time-point scale: {:.8f}".format(scale)
		for p in partials:
			starttime = p[2]
			#print "[{}] starttime: {:.8f}, endtime: {:.8f}".format(p[0], starttime, p[3])
			breakpoints = p[7]
			lasttime = -1
			for bp in breakpoints:
				thistime = bp[0]
				#print "thistime: {:.8f}".format(thistime)
				if thistime > starttime:
					diff = thistime - starttime
					diff *= scale
					bp[0] = starttime + diff
					#print "altering thistime to: {:.8f}, diff: {:.8f}".format(bp[0], diff)
				lasttime = bp[0]
			p[3] = lasttime		# adjust end time
			#print "------------------------------------------"
		return partials

	""" Omit partial if it meets certain criteria, such as freq, amp,
	    duration, and bandwidth ranges. If any breakpoint has an
	    out-of-range value, the entire partial is excluded.
	    breakpoints: tuples of (starttime, freq, amp, [bw])
	    starttime, endtime: spanning all breakpoints for the partial
	    Return None if partial filtered, else return partial breakpoints.
	"""
	def _filter_partial(self, breakpoints, starttime, endtime):
		# exclude partials outside of time range
		if starttime < self._time_range[0] \
				or	(self._time_range[1] != 0 \
					and starttime > self._time_range[1]):
			return None
		if endtime < self._time_range[0]:
			return None
		if	self._time_range[1] != 0:
			if self._clip_end:
				if	endtime > self._time_range[1]:
					numbp = len(breakpoints)
					for i in range(numbp - 1, -1, -1):
						bp = breakpoints[i]
						time = bp[0]
						if time > self._time_range[1]:
							breakpoints.pop()
					if len(breakpoints) < 2:	# can't play these
						return None
					# zero out last bp amp
					breakpoints[-1][2] = 0
			else:
				if	endtime > self._time_range[1]:
					return None

		# exclude partials outside of duration range
		partial_dur = endtime - starttime
		if partial_dur < self._duration_range[0]:
			return None
		if partial_dur > self._duration_range[1] and	self._duration_range[1] != 0:
			return None

		# Loris files can contain zero-duration partials, which we filter.
		if partial_dur <= 0.0:
			return None

		# exclude partials outside of amp range
		(minamp, maxamp) = self.get_partial_amp_range(breakpoints)
		if self._verbose:
			print "minamp: {:.5f}, maxamp: {:.5f}".format(minamp, maxamp)
		if minamp < self._amp_range[0]:
			return None
		if maxamp > self._amp_range[1] and self._amp_range[1] != 0:
			return None

		# exclude partials outside of bandwidth range
		(minbw, maxbw) = self.get_partial_bw_range(breakpoints)
		if self._verbose:
			print "minbw: {:.5f}, maxbw: {:.5f}".format(minbw, maxbw)
		if minbw < self._bw_range[0] or maxbw > self._bw_range[1]:
			return None

		# exclude partials outside of freq range
		(minfreq, maxfreq) = self.get_partial_freq_range(breakpoints)
		if self._verbose:
			print "minfreq: {:.5f}, maxfreq: {:.5f}".format(minfreq, maxfreq)
		if minfreq < self._freq_range[0] or maxfreq > self._freq_range[1]:
			return None

		mean = self.get_partial_freq_mean(breakpoints)

		# exclude partials outside of freq std range
		freqstd = self.get_partial_freq_std(breakpoints, mean)
		if self._verbose:
			print "freqstd: {:.5f}".format(freqstd)
		if freqstd < self._freq_std_range[0]:
			return None
		if freqstd > self._freq_std_range[1] and self._freq_std_range[1] != 0:
			return None

		# exclude partials that are not close enough to a harmonic partial
		# of a user-selected fundamental
		if self._freq_harmonics != None:
			numharms = len(self._freq_harmonics)
			fundamental = self._freq_harmonics[0]
			tolerance = ((self._freq_harm_tolerance * 0.01) * fundamental) / 2.0
			for i in range(numharms):
				p = self._freq_harmonics[i]
				if mean < (p - tolerance):
					return None
				elif mean < (p + tolerance):
					break
				elif i == (numharms - 1):
					return None

		# exclude partials that are not close enough to the pitches of a
		# a user-selected chord
		if self._freq_chord_pitches != None:
			srcpitch = rtcmix.octcps(mean)
			numchordpitches = len(self._freq_chord_pitches)
			mininterval = 999999.0
			index = -1
			for i in range(numchordpitches):
				targetpitch = self._freq_chord_pitches[i] + self._freq_chord_transp
				interval = abs(srcpitch - targetpitch)
				if (interval < mininterval):
					mininterval = interval
					index = i
			if (index < 0 or mininterval > self._freq_chord_sensitivity):
				return None

		# unfiltered partials
		return breakpoints

	""" Helper for _build_freqtab """
	def _make_lfo(self, pdur):
		if self._lfotype != None and pdur >= self._lfodurrange[0] \
				and (self._lfodurrange[1] == 0 or pdur <= self._lfodurrange[1]):
			if self._is_list(self._lforate):
				rate = self._lfo_randgen.choice(self._lforate)
			else:
				rate = self._lforate
			if self._is_list(self._lfomin):
				lmin = self._lfo_randgen.choice(self._lfomin)
			else:
				lmin = self._lfomin
			if self._is_list(self._lfomax):
				lmax = self._lfo_randgen.choice(self._lfomax)
			else:
				lmax = self._lfomax
			if self._lfotype_israndom:
				lfo = rtcmix.makerandom(self._lfotype, rate, lmin, lmax, self._lfoseed)
				if self._lfosmooth > 0.0:
					# Use as initial value for smoothing filter the midpoint between
					# min and max. Otherwise, init value will be zero and will
					# produce a wild glissando when zero is outside of (min,max).
					if self._is_rtcmix_handle(lmin):
						tmin = rtcmix.samptable(lmin, "nointerp", 0)
					else:
						tmin = lmin
					if self._is_rtcmix_handle(lmax):
						tmax = rtcmix.samptable(lmax, "nointerp", 0)
					else:
						tmax = lmax
					initval = tmin + ((tmax - tmin) * 0.5)
					lfo = rtcmix.makefilter(lfo, "smooth", self._lfosmooth, initval)
				self._lfoseed += 1
			else:
				lfo = rtcmix.makeLFO(self._lfotype, rate, lmin, lmax)
			return lfo
		else:
			return None

	""" Helper for _play_one_partial """
	def _build_freqtab(self, envsize, freqpoints, pdur):
		if self._glisstable != None and pdur >= self._glissdurrange[0] \
				and (self._glissdurrange[1] == 0 or pdur <= self._glissdurrange[1]):
			if self._is_list(self._glisstable):
				freqt = self._gliss_randgen.choice(self._glisstable)
			else:
				freqt = self._glisstable
			freqt *= rtcmix.maketable("line", "nonorm", envsize, freqpoints)
		else:
			freqt = rtcmix.maketable("line", "nonorm", envsize, freqpoints)
		lfo = self._make_lfo(pdur)
		if lfo != None:
			freqt *= lfo
		return freqt

	""" Helper for _build_amptab """
	def _make_amplfo(self, pdur):
		if self._amplfotype != None and pdur >= self._amplfodurrange[0] \
				and (self._amplfodurrange[1] == 0 or pdur <= self._amplfodurrange[1]):
			if self._is_list(self._amplforate):
				rate = self._amplfo_randgen.choice(self._amplforate)
			else:
				rate = self._amplforate
			if self._is_list(self._amplfomin):
				lmin = self._amplfo_randgen.choice(self._amplfomin)
			else:
				lmin = self._amplfomin
			if self._is_list(self._amplfomax):
				lmax = self._amplfo_randgen.choice(self._amplfomax)
			else:
				lmax = self._amplfomax
			if self._amplfotype_israndom:
				lfo = rtcmix.makerandom(self._amplfotype, rate, lmin, lmax, self._amplfoseed)
				if self._amplfosmooth > 0.0:
					# Use as initial value for smoothing filter the midpoint between
					# min and max. Otherwise, init value will be zero and will
					# produce a wild change when zero is outside of (min,max).
					if self._is_rtcmix_handle(lmin):
						tmin = rtcmix.samptable(lmin, "nointerp", 0)
					else:
						tmin = lmin
					if self._is_rtcmix_handle(lmax):
						tmax = rtcmix.samptable(lmax, "nointerp", 0)
					else:
						tmax = lmax
					initval = tmin + ((tmax - tmin) * 0.5)
					lfo = rtcmix.makefilter(lfo, "smooth", self._amplfosmooth, initval)
				self._amplfoseed += 1
			else:
				lfo = rtcmix.makeLFO(self._amplfotype, rate, lmin, lmax)
			return lfo
		else:
			return None

	""" Helper for _play_one_partial """
	def _build_amptab(self, breakpoints, envsize, amppoints, pdur):
		if self._ampenv2 == None:
			ampt = rtcmix.maketable("line", "nonorm", envsize, amppoints)
		else:
			if self._is_list(self._ampenv2):
				# atk, atkcrv, hold, relcrv
				atk = self._ampenv2[0]
				atkc = self._ampenv2[1]
				hold = self._ampenv2[2]
				relc = self._ampenv2[3]
				if atk + hold > pdur:
					# try to preserve attack time without a release click
					if atk + 0.01 > pdur:
						atk = pdur - 0.01
						hold = 0
					else:
						hold = min(hold, pdur - atk)
				if hold > 0:
					rel = (pdur - atk) - hold
					tab = rtcmix.maketable("curve", 10000, \
								0,0,atkc, atk,1,0, pdur-rel,1,relc, pdur,0)
				else:
					tab = rtcmix.maketable("curve", 10000, \
								0,0,atkc, atk,1,relc, pdur,0)
				ampt = tab * rtcmix.maketable("line", "nonorm", envsize, amppoints)
			else: # assume _ampenv2 is RTcmix table handle
				ampt = self._ampenv2 * rtcmix.maketable("line", "nonorm", envsize, amppoints)
		amp_static = self._ampscale * 32768.0
		lfo = self._make_amplfo(pdur)
		if lfo != None:
			ampt *= lfo
		return ampt * amp_static

	""" Play one partial.
	    startime, endtime    define partial's original time span
	    breakpoints          breakpoint tuples (time, freq, amp) for the partial
	    Return a tuplet of partial start and end times, after timescaling --
	    in other words, the params passed to the RTcmix instrument.
	"""
	def _play_one_partial(self, pindex, numpoints, starttime, endtime, phase, breakpoints):
		start = breakpoints[0][0]
		if start != starttime:	# does this ever happen?
			print "WARNING: First breakpoint time is not the same as partial start time (for partial {}).".format(pindex)
		end = breakpoints[-1][0]
		if end != endtime:
			print "WARNING: Last breakpoint time is not the same as partial end time (for partial {}, endtime: {:.8f}, end: {:.8f}).".format(pindex, endtime, end)
		#framedur = breakpoints[1][0] - start	# not valid if time reassignment performed during analysis
		dur = endtime - starttime
		if dur <= 0.0:	# should already been trimmed during _filter_partial
			return
		tscale = 1.0
		if self._is_number(self._timescale):
			tscale = self._timescale
		else:
			# FIXME: This is not even remotely correct. To do this right, would
			# have to copy a section of the table between starttime and endtime
			# and store *that* in a table. Then use this table to warp the timing
			# of the freq and amp tables, as well as warp the start time and
			# duration of the note. 
			print "calling rtcmix.tablelen; tscale type is", type(self._timescale)
			tablen = rtcmix.tablelen(self._timescale)
			index = (start / self._totdur) * tablen
			tscale = rtcmix.samptable(self._timescale, index)
			#print "index:", index, "tscale:", tscale
		dur *= tscale
		numpoints = len(breakpoints)
		if self._verbose:
			print "[{}] -----------------------------------------".format(pindex)
			print " numpoints: {}, starttime: {:.5f}, endtime: {:.5f}, dur: {:.5f} (post time-scale)".format(numpoints, starttime, endtime, dur)
		if numpoints > self._max_num_breakpoints: 
			print "WARNING: Partial #{:d}, beginning at {:.5f} seconds with {:.2f} Hz, requires more breakpoints ({:d}) than RTcmix can handle ({:d}). Skipping it.".format(pindex, start, breakpoints[0][1], numpoints, self._max_num_breakpoints)
			return (0, 0)

		# build breakpoint time,value temp arrays
		freqpoints = []
		amppoints = []
		bwpoints = []
		bw_has_nonzero = False
		for bp in breakpoints:
			time = bp[0] - start
			freq = bp[1]
			if freq < 15.0:
				freq = 15	# otherwise interpreted as oct.pc by WAVETABLE
			amp = bp[2]
			bw = bp[3]
			if bw > 0:
				bw_has_nonzero = True
			#print "time:", time, "freq:", freq, "amp:", amp, "bw:", bw
			freqpoints.append(time)
			freqpoints.append(freq)
			amppoints.append(time)
			amppoints.append(amp)
			bwpoints.append(time)
			bwpoints.append(bw)
		if self._zero_first_partial_amp:
			amppoints[1] = 0.0
		if self._zero_last_partial_amp:
			amppoints[-1] = 0.0
		if self._zero_phase:
			phase = 0
		if self._verbose:
			print "starttime:", starttime, "endtime:", endtime, "phase:", phase
			print "freqpoints:", ['{:.5f}'.format(item) for item in freqpoints]
			print " amppoints:", ['{:.5f}'.format(item) for item in amppoints]
			print "  bwpoints:", ['{:.5f}'.format(item) for item in bwpoints]
		envsize = numpoints * self._kEnvSlotsPerPoint

		#print "making freq env"
		freqt = self._build_freqtab(envsize, freqpoints, dur)
		#print "making amp env"
		ampt = self._build_amptab(breakpoints, envsize, amppoints, dur)
		if bw_has_nonzero:
			#print "making bw env"
			bwt = rtcmix.maketable("line", "nonorm", envsize, bwpoints)
		else:
			bwt = 0
		if self._timescale_start_points:
			start *= tscale
		pan = 0
		if self._nchans == 2:
			pan = self._pan_randgen.uniform(0, 1)

		# play note on either BWESINE (if necc.) or WAVETABLE
		if bw_has_nonzero or phase != 0.0:
			if self._nchans > 2:
				ch = self._pan_randgen.randint(0, self._nchans - 1)
				if self._multichan_bustype is "out":
					rtcmix.bus_config("BWESINE", "out " + str(ch))
				else:
					rtcmix.bus_config("BWESINE", "aux " + str(ch) + " out")
			rtcmix.BWESINE(start, dur, ampt, freqt, bwt, phase, pan, self._wavetable)
		else:
			if self._nchans > 2:
				ch = self._pan_randgen.randint(0, self._nchans - 1)
				if self._multichan_bustype is "out":
					rtcmix.bus_config("WAVETABLE", "out " + str(ch))
				else:
					rtcmix.bus_config("WAVETABLE", "aux " + str(ch) + " out")
			rtcmix.WAVETABLE(start, dur, ampt, freqt, pan, self._wavetable)
		return (start, start + dur)

	""" Play all unfiltered, unmuted partials.
	       partials     list of partials to play [default is all]
	    Return a tuplet of earliest start time and latest end time.
	"""
	def play_partials(self, partials=None):
		if partials == None:
			partials = self._partials
		minstart = sys.float_info.max
		maxend = 0
		for partial in partials:
			mute = partial[6]
			if not mute:
				pindex = partial[0]
				numpoints = partial[1]
				starttime = partial[2]
				endtime = partial[3]
				phase = partial[4]
				breakpoints = partial[7]
				(start, end) = self._play_one_partial(pindex, numpoints, starttime, endtime, phase, breakpoints)
				if start < minstart:
					minstart = start
				if end > maxend:
					maxend = end
		return (minstart, maxend)

# utilities ###################################################################

""" Given a list or tuple containing pitch classes in oct.pc format, return
	 a list of pitches in specific registers from minoctave up to, but not
	 including, maxoctave. This lets you specify a scale or collection and then
	 arpeggiate it across multiple octaves. For example,

	    pc = (0.00, 0.02, 0.04, 0.05, 0.07, 0.09, 0.11)   # major scale
	    pitches = pclist2pitchlist(pc, 8, 10)
"""
def pclist2pitchlist(pclist, minoctave=4, maxoctave=14):
	pitchlist = []
	if not (minoctave < maxoctave):
		print "\nERROR: pclist2pitchlist: minoctave must be less than maxoctave."
		return None
	for pc in pclist:
		for octave in range(minoctave, maxoctave):
			pitchlist.append(octave + pc)
	return sorted(pitchlist)


# for testing #################################################################
if __name__ == "__main__":
	sr = 44100
	nchans = 2
	writeit = 0
	#infile = "_test/boysong-excerpt.txt"
	infile = "/Users/johgibso/pieces/edgewater/par/531-04-R.txt"
	totdur = 15.0
	outfile = "test.wav"
	if writeit:
		rtcmix.set_option("clobber=yes", "play=no")
		rtcmix.rtsetparams(sr, nchans)
		rtcmix.rtoutput(outfile, "24")
	else:
		rtcmix.rtsetparams(sr, nchans)
	#rtcmix.control_rate(sr)
	rtcmix.control_rate(sr / 4)

	rtcmix.set_option("print = 1")
	r = Resynth(sr, nchans)
	#r.select_partial_range(0, 1000)
	r.select_time_range(1, 4)
	#r.select_duration_range(0.002, 0.02)
	#r.select_amp_range(0.0001, 0.02)
	#r.select_freq_range(0, 0)
	#r.select_freq_max_std(1.4)
	fund = rtcmix.cpspch(8.00)
	#r.select_harmonic_partials(fund, 20)
	#r.select_bandwidth_range(0, 0)
	#r.set_verbose(True)
	numpartials = r.read_file(infile)
	if numpartials > 0:
		r.set_zero_first_partial_amp(True)
		r.set_zero_last_partial_amp(True)
		#r.set_gain(0.0)
		r.clamp_bandwidths(0.0, 0.2)
		#r.scale_bandwidths(0.5)
		r.scale_bandwidths(rtcmix.maketable("curve", "nonorm", 1000, 0,0,5, 1,2))
		#r.scale_and_offset_freqs(0.5, 30, 0.002)
		chord = (7.00, 7.07, 8.00, 8.02, 8.07, 9.00, 9.02, 9.03, 9.07, 10.00)
		chord = (6.10, 7.00, 7.04, 7.05, 7.07, 7.10, 8.00, 8.04, 8.05, 8.07, 8.10, 9.00, 9.02, 9.05, 9.07, 10.00)
		#r.retune_partials(chord, transp=2, sensitivity=36.0, strength=0.9)
		#r.synthesize_new_partials(100, sr/4, 3.0, -9.0, 0.002)
		#r.set_wavetable(rtcmix.maketable("wave", 32767, "tri5"))

		(starttime, endtime) = r.get_time_bounds()
		print "time bounds: ({:.5f}, {:.5f})".format(starttime, endtime)
		timescale = 100
		#r.scale_time_points(0.005)
		#r.set_timescale(timescale)
		#r.set_timescale_start_points(False)
		timeshift = -starttime
		r.shift_times(timeshift)

		r.set_pan_seed(3)
		r.set_outbus("aux 0-1 out")
		rtcmix.bus_config("FREEVERB", "aux 0-1 in", "out 0-1")
		rtcmix.load("FREEVERB")
		mix = 10
		rtcmix.FREEVERB(0, 0, (totdur + timeshift) * timescale, 1, 0.7, 0.1, 2, 0, 100-mix, mix, 100)
		r.play_partials()
		#r.write_file("boyout.txt")
