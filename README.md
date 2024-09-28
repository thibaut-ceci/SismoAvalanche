# What is "SismoAvalanche"?

SismoAvalanche is a project with two objectives :
- Try to find the morphometric parameters of avalanches using the seismic signal
- Improve the understanding of avalanche collapse dynamics by exploring a new field : avalanche entropy

The report (in french) of this study is available here and it is advisable to read it before consulting the codes : https://mega.nz/fm/cNVHSLoa

We will use the database ESEC : https://ds.iris.edu/spud/esec

## Abstract of this study :
Understanding the seismic characteristics of landslides is crucial for assessing their impact and mitigating risks. This study focuses on using seismic signals to deduce morphometric parameters (volume, height of fall and runout) of dense granular flows, such as debris, snow, ice, or rock avalanches. Seismic data from 30 cataloged avalanches in the Exotic Seismic Events Catalog (ESEC) are downloaded and analyzed to extract and compare twelve characteristics from the power spectrum, energy envelope and entropy spectrum of these events. The results reveal interesting links, proportional relationships and correlations : the morphometric parameters of avalanches show proportional relationships and the characteristics of the energy envelope allow for an approximate estimation of the volume of these events. However, the characteristics based on entropy did not show significant correlations or links with other characteristics. This is partly caused by several issues related to the data provided in the ESEC and the detection method that leads to the loss of many events constraining the study to a reduced sample of avalanches. With more events, links or correlations could be observed. Furthering this study could improve the understanding of dense granular flows.


# Explanation of the codes :
- 00_read_esec.ipynb : First, we will read the ESEC and only keep the avalanches (because in ESEC, there are a lot of types of landslides, not only avalanches).
- 01_clean_esec.ipynb : Then, we will remove unnecessary columns in the new ESEC and keep only events with good measurements.
- 02_compute_incertainties.ipynb : Not all events have uncertainties. We will fill them.
- 03_download_seismic_waveform.ipynb : Seismic waveforms of the kept events will be downloaded.
- 04_clean_pickle.ipynb : Unnecessary pickle files will be removed (like pickle files without instrumental response).
- 05_explore_esec.ipynb : Here, we will explore the ESEC and compute the number of stations for each event.
- 06_event_mapper.ipynb : In this notebook, we want to see the location of the avalanches.
- 07_detection_one_trace.ipynb : A detection method was created in this study. It is tested in this notebook.
- 08_detection_all_events_stream.ipynb : The detection method was applied to all events. A plot was created for each event.
- 09_detection_all_events.ipynb : The detection method was applied to all events again. A plot was created for all events.
- 10_fitting_spectrum.ipynb : Using the Welch method, a spectrum was created for each event. To extract the features, two models were created and fitted to the spectrum.
- 11_compute_energy.ipynb : Using the spectrogram, the energy envelope was computed and some features were extracted.
- 12_compute_entropy.ipynb : The entropy of each event was computed and some features were extracted.
- 13_features_comparison : All the features computed in this study were merged into a dataframe.
- 14_scatter_matrix : A scatter matrix was created to look for correlations or links between the features.

# A lot of librairies developped in this studies were availables :
- analysis.py : Analyses seismic waveforms with ObsPy. The detection method is here.
- catalog.py : ESEC catalog management.
- cleaning.py : Libraries to clean ESEC and associated datas (pickle files).
- computations.py : ESEC catalog computation libraries.
- energy.py : Performs energy calculations.
- entropy.py : Performs entropy calculations.
- figures.py : Manages the creation of figures.
- mapper.py : Manages the creation of maps.
- matrix.py : Create a scatter matrix.
- waveform.py : Downloads the seismic data.

# Explanation of the files :
- \catalog : Contains the many catalogs created during this study
- \features : Contains the features calculated during this study
- \figures : Contains the plots obtained during this study


This study took place between February 2024 and October 2024 as part of Thibaut Céci's Master 2 internship, which was supervised by Antoine Lucas (lucas@ipgp.fr), Léonard Seydoux (seydoux@ipgp.fr) and Anne Mangeney (mangeney@ipgp.fr).
