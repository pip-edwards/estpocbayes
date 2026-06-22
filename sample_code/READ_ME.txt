----------------------------
input_data_setup.py

	set up the SST, Chl-a and depth into matched spatial grids
	make climatologies for SST and Chl-a	

----------------------------
match_poc.py

	matches the poc database to satellite derived sst and chl-a data
	makes figure S1 for showing map of data availability
	split into 10%test 90% train data

----------------------------
stan_run.py

	runs the stan code (stan_noy.stan) on the existing database
	exports a table of the summarised stats
	saves all outputs including beta gamma mu and sigma separately into its own table

---------------------------

histo_stats.py

	plots Figure 2 and creates a table of comaprison statistics used in supplementary

------------------------------

effect_sizes.py

	plots Figure 3

------------------------------

global_poc.py

	plots maps of poc flux and attenuation rate (Figure 4), (S6-S8)
	makes global estimate of poc flux. 
	WARNING last bit is very very unoptimized loop. I wrote to run overnight/while writing and did not bother to improve :D

----------------------------
