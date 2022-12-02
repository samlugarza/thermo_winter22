from astropy.table import Table
dat = Table.read('/Users/samgarza/Documents/thermo_winter22/pset2/firas_monopole_spec_v1.txt', format='txt')
print(dat)