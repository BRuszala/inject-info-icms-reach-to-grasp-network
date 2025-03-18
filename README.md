# inject-info-icms-reach-to-grasp-network

The goal of this project was to determine in which of 5 cortical regions in the cortical reach-to-grasp network information could be injected via low-amplitude ICMS: primary somatosensory cortex (S1), ventral premotor cortex (PMv), dorsal premotor cortex (PMd), anterior intraparietal area (AIP), or dorsal posterior pareital cortex (dPPC).  

This repository contains .h5 fiels with all task parameters and behavioral data, .ns2 files with analog signals of the joystick movements and curosor position of the screen sampled at 1 kHz, and .nev files that contain event markers to identifiy individual epochs of the task:

Event Marker Values:
30	prior trial was an error	
31	prior trial was a success	
60	reach to target	
85	Initial wait	
100	error	1 byte follows: error description code
119	Success	
120	Reward(s) on	
121	Reward 1 off	
122	Reward 2 off	
210	instruction cue	
211	delay	
245	go cue	
248	target hold start	
250	pre inter trial	
251	Trial start	12 bytes follow; see table below
253	Trial end

"Trial Kind" information: 
0	  all errors enabled
1	  wrong target error disabled
2	  reach to center error disabled
3	  = 1 + 2
4	  = final hold error disabled
5	  = 1 + 4
6	  = 2 + 4
7	  = 1 + 2 + 4
10	= 0 + catch trial
11	= 1 + catch trial
12	= 2 + catch trial
13	= 3 + catch trial
14	= 4 + catch trial
15	= 5 + catch trial
16	= 6 + catch trial
17	= 7 + catch trial

See an example trial:
1	  251	start of trial
2	  serno_1	231	serial number byte 1; serial_number = serno_1 + serno_2*256 + serno_3*256^2 + serno4*256^3 = 609903
3	  serno_2	58	serial number byte 2
4	  serno_3	9	serial number byte 3
5	  serno_4	0	serial number byte 4
6	  255	unused
7	  trial_kind	0	see above; reach to center error disabled
8	  annuli_count	1	annuli count in current task
9	  targets_count	8	targets count in current annulus
10	selected_annulus	0	selected annulus index (0-based)
11	selected_target	2	selected target index (0-based, counterclockwise where target in +X direction is 0)
12	prior_success_status	31	30 (prior trial was an error) or 31 (prior trial was a success)
13	trial_execution_status	34	34 (trial is an execution trial) or 35 (trial is an observation trial)
14	210	instruction cue
15	245	go cue
16	60	60	reach (left center)
17	248	entered target; target hold start
18	119	success (target hold end)
19	120	reward on
20	121	reward off
21	250	pre inter trial; just extra period for solenoid relaxation)
22	253	end of trial

EMG Data are also included for constructing stimulus-triggered averages.
