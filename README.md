# signcor

This repository contains the code and the data used in the paper "Estimating correlations among elliptically distributed random variables under any form of heteroskedasticity".

The files are organized as follows.

-   **data** folder
    -   `SP100_constituents_5min.Rdata` 5-minute prices of SP100 consituents.
    -   `SP100_30y.xlsx` daily prices of SP100 constituents.
-   **code** folder
    -   `application.R` Code to generate the plots contained in Figure 5. To change the stock pairs you have to reassign the string variables `stck1` and `stck2`. It uses the R data file `SP100_constituents_5min.Rdata`.
    -   `dyn_cor_sim.R` Code to carry out the simulations and application of Section 4. It uses the data in the Excel file `SP100_30y.xlsx`.
    -   `plots_maker.R` Code to produce Figure 9 and 10.
    -   `sim_rho_estim.R` Code to produce the simultions of Section 2.2. Version with known medians.
    -   `sim_rho_estim_demean.R` Code to produce the simultions of Section 2.2. Version with estimated medians.
