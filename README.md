## Python and Julia code to support manuscript:
### D. Heslop, U. Amarathunga and E. J. Rohling (2023). North African Plio-Pleistocene monsoon runoff and its impacts on the Mediterranean Sea.

Code and data files are provided in two folders to recreate manuscript Figures 2, 3, 5, 6, 7 and 8.

### Quantile regression
A Jupyter notebook ```Quantile regression - ODP967.ipynb``` with Python code is provided to recreate Figures 3, 4, 5, and 6. The notebook contains code comments and a detailed description of the processing steps. The folder also contains the following data sets needed for plotting:
1. ```ODP967 - d18O.csv``` - Tuned ODP Site 967 $\delta^{18}\mathrm{O}_{\mathrm{pf}}$ record.
2. ```ODP967 - Nonsapropel ages.csv``` - Ages of nonsapropel intervals in the ODP Site 967 record.
3. ```Rohling2021 - ESL.csv``` - Eustatic sea-level record of Rohling *et al.* (2021, doi:10.1126/sciadv.abf5326)
4. ```Laskar2004 - Eccentricity.csv``` - Eccentricity record of Laskar et al. (2004, doi:10.1051/0004-6361:20041335)

When running code locally, you will need to include path information in the ```np.loadtxt``` statements if the listed data files are not in the same folder as the notebook.

### Box model sensitivities
Contains Jupyter notebook ```Box model sensitivities.ipynb``` with Julia code to estimate the Winter-Summer and Winter-Summer-Monsoon box model sensitivities in Figures S2 and S4.  The notebook contains code comments and a detailed description of the processing steps. The provided Julia code files:
1. ```rsl.constant.jl```
2. ```rsl.estimate.jl```
3. ```rsl.fixvar.jl```
4. ```rsl.forcing.jl```
5. ```rsl.rndvar.jl```
6. ```rsl.utils.jl```

are also required and should be placed in the same folder as the notebook.

Once the calculations in ```Box model sensitivities.ipynb``` are complete, the resulting output files can be used to plot Figures 2 and 8 using the Python code in ```Box model sensitivities plotter.ipynb```. You should ensure that the ```Box model sensitivities.ipynb``` output files are in the same folder as the ```Box model sensitivities plotter.ipynb``` notebook. If not, you will need to include path information in the ```np.loadtxt``` statements in ```Box model sensitivities plotter.ipynb```. 

 
