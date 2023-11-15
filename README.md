# ANdataset
This repository contains matlab scripts that are used to navigate and analyzed data from auditory nerve fiber recordings. 

**Software to search through the dataset**

Three scripts are provided to help the user search through the full dataset as well as within a struct of one animal. 
1.	check_dataset_metadata.m loops through the ‘all_exp’ struct and focusses on the metadata of the experiments. It recreates the metadata sheet, which can be used to select an animal of interest and investigate it further in check_dataset_animal.m. This script was used to generate Table 1 and Figure 2.
2.	check_dataset_units.m loops through the ‘all_exp’ struct and focusses on the analysed outcomes of the single units. It can be used to plot any of these outcomes against each other, typically with the fibre’s BF on the x-axis. check_dataset_units.m was used to generate Figure 6a (BF vs threshold), Figure 6b (BF distribution), Figure 6c (BF vs SR), Figure 6d (BF vs VS), and Figure 5b (BF vs click latency).
3.	check_dataset_animal.m loops through the units of an ‘exp’ struct of one animal. It generates a scatterplot of the threshold against BF of all the fibre’s recorded in that animal and a plot with all the RLFs of the animal plotted in one graph (e.g. Figure 5a, generated from ‘G220922.mat’). This script calls the function check_AN.m, which plots the first trial of the recording, the inter-spike interval histogram, the first 300 spike waveforms, and the median spike waveform +/- 95% confidence interval. The output of this function is shown in Figure 4 (BF recording of ‘G220922.mat’, unit ‘3p_607’ [i = 12]). The input is a curvedata, curveresp, and curvesettings field from the same recording. It also calls the function makePSTH.m, which is used to generate a peri-stimulus time histogram of all responses to tone bursts at or close to a given stimulus level above the fibre’s threshold. makePSTH.m was used to generate Figure 5c, based on the recorded spike times of animalID ‘G220908’ from unit ‘3p_181’ (i = 23) at 20 dB above threshold (TestLevel = 20).

**Matlab code that was used to generate the analysed outcomes in the data structs**

The script call_extract_func.m calls all functions that were used to generate the analysed outcomes of one unit of one animal. These include the following functions:
1.	BFextract_func.m, which generates the frequency-response curve to derive the BF.
2.	CFextract_func.m, which generates the receptive field and tuning curve to derive the characteristic frequency, threshold, and Q10dB. 
3.	PHextract_func.m, which calculates the vector strength and plots it as a function of stimulus level. 
4.	CLICKextract_func.m, which calculates the click latency in three different ways and visually demonstrates the outcomes of each method. 
5.	RLFextract_func.m, which generates the rate-level function to derive the threshold. 
6.	SRextract_func.m, which calculates the spontaneous rate from a long recording in silence. 

