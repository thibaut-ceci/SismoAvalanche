
A mettre à jour



# Analyse seismic signals of avalanches

This is the repository of seismic signal analysis for debris/rock/snow avalanches.

## Getting started

We gather event catalog from https://ds.iris.edu/spud/esec

Then, we select event which correspond to avalanche (i.e., catalog = catalog[catalog["subtype"].str.contains("avalanche")])

The `get_waveforms_from_catalog.ipynb` provides some starting point to get waveforms for each event in the catalog. The may take a while to process accross the whole data set.

The `get_entropy.py` is a starting notebook that simply applies the coherence of the covariance over the trace in order to derive entropy of the seismic signal. The expected result is that antropy goes down during an avalanches event. This is shown in the example. 


## On the analysis of the seismic signals

- [x] Get waveforms from catalog
- [x] Turn into functions the cell of the spectrogram
- [ ] Calculate the source time function parameters from the source spectra