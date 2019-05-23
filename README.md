# grangers
This repository contains data and code relative to the manuscript "A bootstrap test to detect prominent Granger-causalities across frequencies" by Matteo Farnè and Angela Montanari (https://arxiv.org/abs/1803.00374).

We have created an R package, called "grangers", containing seven functions:
- "Granger.unconditional", performing the calculation of the unconditional Granger-causality spectrum;
- "Granger.conditional", performing the calculation of the conditional Granger-causality spectrum;
- "Granger.inference.unconditional", performing boostrap inference on the unconditional Granger-causality spectrum;
- "Granger.inference.conditional", performing boostrap inference on the conditional Granger-causality spectrum;
- "Granger.inference.difference", performing boostrap inference on the difference between an unconditional 
   and a conditional Granger-causality spectrum;
- "bc_test_uncond", performing the parametric test of Breitung and Candelon (2006) on the unconditional Granger-causality spectrum;
- "bc_test_cond", performing the parametric test of Breitung and Candelon (2006) on the conditional Granger-causality spectrum.
   
For the details we refer to the help of each function and to the paper (https://arxiv.org/abs/1803.00374). 
   
The R package "grangers" also stores the dataset "data_all", which contains the time series used in the analysis (GDP, M3, M1, HICP, UN, LTN) for the periods 1999:Q1-2017:Q4 respectively. For the details we refer to the paper (https://arxiv.org/abs/1803.00374) and to  http://sdw.ecb.europa.eu/web/generator/docu/rtdb_docu.pdf.
