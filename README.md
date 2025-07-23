# Beyond sex differences in mean: meta-analysis of differences in skewness, kurtosis, and covariance

html file contains the analyses necessary to reproduce this research. Alternatively, one can check the rmd file for the code. Data comes from a single csv file:

## `mice_data_sample.csv`

Sample from the International Mouse Phenotyping Consortium (IMPC; Dickinson et al., 2016). Each row contains information for a different mouse for a specific trait. 

### metadata
- `specimen_id`: unique identifier for articles in our dataset; open-ended text
- `strain_accession_id`: mice strain; open-ended text
- `phenotyping_center`: phenotyping center; open-ended text
- `value`: trait value; numeric
- `trait_name`: name of the trait measured; open-ended text

### Notes
- Formulas for analytical approximation for sampling variances for skewness and kurtosis are seem to only apply to situations where data are normally distributed.