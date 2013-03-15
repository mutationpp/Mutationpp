clear

nCH4=.067;
nCO=.037;
nCO2=.011;
nH2=.334;
nH2O=.489;
nC6H5OH=.062;
nCH4+nCO+nCO2+nH2+nH2O+nC6H5OH

nC=nCH4+nCO+nCO2+6*nC6H5OH;
nH=4*nCH4+2*nH2+2*nH2O+6*nC6H5OH;
nO=nCO+2*nCO2+nH2O+nC6H5OH;
sum=nC+nH+nO
nC=nC/sum
nH=nH/sum
nO=nO/sum
