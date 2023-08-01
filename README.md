# redneb

Calculate nebular extinction from hydrogen Balmer series emission lines, 
accounting for temperature and density of the gas. Stellar absorption optionally 
included in the solution.

## Example calling sequence

```
from redneb import *
# example Balmer wavelengths
w = np.array([6563,4861,4341,4102])
# example fluxes
f = np.array([471.87,145.13,66.31,35.22])
# example flux errors
e = np.array([6.67,3.06,1.81,1.56])
# example equivalent widths -- used only for stellar absorption
ew = np.array([978.20,160.47,72.64,37.91])
# keyword arguments for reddening functions, assuming temp and density
kwarg = dict(ext_law='ccm89',flux_err=e,temp=1e4,dens=1e2)
# E(B-V) with no stellar absorption
ebv_arr,ebv_err = calc_ebmv_from_balmer(w,f,**kwarg)
ebv_avg = np.average(ebv_arr,weights=1/ebv_err**2)
# E(B-V) with stellar absorption
ebv_fit = fit_ebmv_starabs(w,f,ew,**kwarg)[0]
# print results
print(f'E(B-V) without star absn: {ebv_avg:.3f}')
print(f'E(B-V) with star absn: {ebv_fit:.3f}')
```

## Licensing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
