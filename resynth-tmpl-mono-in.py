import rtcmix, resynth
#maxdispargs = 4096
rtcmix.set_option("print = 1")

parfile_name = "lake-waves.txt"
outfile_name = "resynth-tmpl-mono-in.wav"
#paroutfile_name = "resynth-tmpl-mono-in-par.txt"

sr = 44100
cntlrate = sr
numoutchans = 2   # from Resynth to out or to stereo reverb
use_aux_out = 0   # only for multichan out
writeit = 0
usereverb = 1     # not for production output
reverb_mix = 5    # wet mix %
reverb_roomsize = 0.7

verbosity = 0     # 0: quiet, 1: resynth.py, 2: rtcmix, 3: resynth.py + rtcmix
# print params for all selected partials within time and freq ranges; don't play
histogram = 0     # 0: no, 1: print before processing
printparams = 0   # 0: no, 1: print before processing, 2: print after processing

gain = 24

### selections (comment out to ignore selection) ############################
selstart = 0; selend = 0
minfreq = 60; maxfreq = 0
mindur = 0.1; maxdur = mindur + 0.15
#minamp = 0.0001; maxamp = 0.2
#minbw = 0.1; maxbw = 0.9
# select partials whose frequency points deviate by no more than std...
#minstd = 2; maxstd = 15
# select harmonic partials above a given fundamental...
#fundamental = rtcmix.cpspch(8.00); tolerance = 20  # tolerance is % of fund.

### modifications (comment out to skip) #####################################
seed = 1

#-- time --------------------------------------
timescale = 2
#timescale_start_points = False
#additional_timepointscale = 1   # useful when timescale_start_points is true
#delaytab = rtcmix.maketable("line", "nonorm", 10000, \
#	0,0, 1000,2, sr/2,4); deldev = 0.0
#quantum = 0.25; quantdev = 0.00; quantsmear = 0.01; quantseed = 3; quanttmpl = (1, 1, 1, 0, 1, 0, 1)

#-- amplitude ---------------------------------
#clampminamp = 0.01; clampmaxamp = 0.01
#invert_amps = True
#zero_first_amp = True   # within each partial
#zero_last_amp = True    # within each partial
#envtab = rtcmix.maketable("curve", 1000, 0,1,1, 0.01,1,-4, 1,0)
envtab = (0.0005, 0, 0, -8)

#-- frequency ---------------------------------
#retune_last = True
#freqscale = 1.0; freqshift = 0; freqdev = 0.00
#freqscale = rtcmix.maketable("curve", "nonorm", 10000, 0,0.3,-1, 1,1,0, 2,1); freqshift = 0; freqdev = 0
#chord = pclist2pitchlist((0.00, 0.02, 0.03, 0.05, 0.07, 0.10), 6, 14)
transp = 0; sensitivity = 4.0; strength = 1.0
#strength = rtcmix.maketable("curve", 1000, 0,1,4, 5,0,0, 6,0)
# synthesize new partials...
#sminfreq = sr/9; smaxfreq = sr/4; smult = 2; sgain = -3; sdev = 0.005
#wavet = rtcmix.maketable("wave", 32767, "tri5")
#eqtab = rtcmix.maketable("line", "nonorm", 10000, \
#	0,0, 200,-12, 800,-12, 1600,6, 4000,0, sr/2,-18)
# make a bunch of tables that will be randomly chosen to shape freq curves
#glisstab = []
if 'glisstab' in globals():
	import random
	random.seed(999)
	gdur = (0.5, 4.0)   # range of partial durations to affect
	for i in range(0, 20):
		a = random.uniform(0.8, 1.2)
		b = random.uniform(0.8, 1.2)
		t = rtcmix.maketable("curve", "nonorm", 1000, 0,1,3, 1,a,3, 2,b)
		glisstab.append(t)
#lfotype = "square"
if 'lfotype' in globals():
   lrate = rtcmix.maketable("line", "nonorm", 1000, 0,4, 1,12, 2,2)
   lmin = rtcmix.maketable("line", "nonorm", 1000, 0,1, 1,0.95, 2,1)
   lmax = rtcmix.maketable("line", "nonorm", 1000, 0,1, 1,1.05, 2,1)
   lseed = 1
   ldmin = 0.5; ldmax = 0

#-- bandwidth ---------------------------------
#clampminbw = 0.1; clampmaxbw = 1.0
#scalebw = 0.0
#bwcurve = rtcmix.maketable("curve", "nonorm", 1000, 0,0,0, 0.7,0,2, 2,1)
#bwcurve = rtcmix.maketable("linestep", 1000, 0,0, 1,1, 2,1)

#==============================================================================
# execute common code using params defined above
f = open("./resynth-common-mono.py")
code = f.read()
exec(code)
f.close()

