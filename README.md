Code for the paper:

Federated Ensemble Regression using Classification. DS 2020.

To reproduce the results having copied to repository to a local directory:
1. Download the [datasets](http://dx.doi.org/10.17632/8mgyb6dyxv.2) and place them in /input/datasets/
2. Execute feruc.R in the code directory
3. The performance for both the regressors and classifiers are logged in /output/performance/. Models are logged in /output/models/ 

_Note_: The output performance files are named as follows 'experiment type'--'number of bins'--'aggregator'--'discretizer'.txt. Experiment type can either be "cls" for classification or "reg" for regression. Number of bins is an integer within the limit specified in the paper (2-5). Aggregator can be "ag" for simple averaging, "dn" for undersampling, "up" for oversampling, and "rg" for the case in which class imbalance is ignored. Discretizer can be "freq" for equal frequency, "kmn" for k-means, "srt" for simple sorting, and "rdm" for random binning.

In the regression performance files, the first row is the name of the gene, followed by the R squared, mean square error, and root mean square error. 
The classification performance found in /output/performance/ is the mean performance for each gene which can be found in /output/performance/all_class. These results follow the above naming convention, however, "cls" is replaced with the gene name.

Given the number and size of models, it might be useful to disable model logging. This can be done by commenting out line 284 (`log_models(models, gene, nbins, type)`).

The experiments were performed in R 3.4.4, and require the following libraries: data.table, ranger, smotefamily, caret and arules.
