# README

Related [GitHub Issue](https://github.com/scverse/mudata/issues/111)

## Description

> Build a scverse-native data format that accounts for the hierarchical nature of quantification in MS-based proteomics.

### Problem Statement

In LC/MS-proteomics, there is a naturally arising hierarchical feature structure

- At the lowest level, mass spectrometers detect + quantify _fragments_ from charged peptides (_precursors_) in the mass spectrometry (MS) instruments. The precursor-level data is relatively large (N samples x ~100 000 features)
- Proteomics search engines identify the precursor sequences and match them to their corresponding proteins. Ultimately, search engines derive _protein_-specific intensities (N samples x 3000-10000 features).

The key challenge is that there exists an N:M relationship between precursors and proteins, i.e. many precursors can map to one protein and sometimes a precursor could potentially be derived from different (homologous) proteins.

The main extension of the data format to existing data containers like mudata would be the formalization the relationship/mapping between the fundamental units of quantification in MS-proteomics (fragments, precursors) and high-level, biologically more relevant aggregated features (peptides, proteins, genes).

- [ ] Implement an RFC (see [rfc/RFC.md](rfc/RFC.md)). 
- [x] Implement the prototypes of the data structure that have been proposed in an scverse (i.e. anndata/mudata) compatible manner. See [`src/`](src/) and [`notebooks/demo_hierarchical_mudata.ipynb`](notebooks/demo_hierarchical_mudata.ipynb).
- [x] **Application**: Implement a related, simple downstream analysis that builds on the data format to get an intuition for the API (e.g. “Plot the distribution of all precursor intensities that correspond to a specific protein”)
- [ ] **Data ingestion**: Implement one proof-of-principle reader from a quantification pipeline/search engine output (e.g. QuantMS, DIANN, alphadia) to the data container.

## Support 

## Get started

We recommend contributors to make themselves familiar with the [mudata](https://mudata.readthedocs.io/stable/notebooks/nuances.html) documentation and API. 

See also [QFeatures](https://rformassspectrometry.github.io/QFeatures/articles/QFeatures.html) for a conceptually similar R-package and [alphaquant](https://github.com/MannLabs/alphaquant.git) for a potential future application of the mapping approach.

## Example data 
Use the download script to obtain example PSM reports. 

```shell
bash download.sh

# download real world data
bash data/download.sh -o data/ albrecht2025

# download a minimal dataset 
bash data/download.sh -o data/ minimal
```

