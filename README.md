# Description
This repository contains the fortran scripts allowing to generate the figures of the paper "Control analysis of cooperativity and complementarity in metabolic regulation; the case of NADPH homeostasis" by Pfeuty et al. , metabolites, 2023.

The repository includes 10 files:
- fig2.f90 and fig2.out are the fortran code and the executable to generate data plotted in Figure2. The code fig2.f90 reads the file param*.dat that contains the parameters associated to various subdomains of the flux space.
- fig3.f90 and fig3.out are the fortran code and the executable to generate data plotted in Figure3.
- fig4.f90 and fig4.out are the fortran code and the executable to generate data plotted in Figure4.

The compilation requires three fortran libraries : lapack minpack and seulex.
