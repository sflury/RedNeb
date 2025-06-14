# RedNeb

Calculate nebular extinction from hydrogen Balmer series emission lines, 
accounting for temperature and density of the gas. Stellar absorption 
corrections are optionally included in the solution.

While this code is provided publicly, I request that any use thereof be 
cited in any publications in which this code is used. 
The publication which originally introduced and used this code is
Flury et al. 2022 ApJS 260, 1. Please also cite any relevant
relevant extinction laws (e.g., Cardelli et al. 1989. ApJ 345, 245),
the PyNeb package (Luridiana et al. 2015. A&A 573, 42), and, if applicable,
the stellar absorption formalism (Izotov et al. 1994. ApJ 435, 647)

## Example usage

``` python
from RedNeb import *
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

## BibTeX reference (as provided by NASA/ADS)
Flury et al. 2022 ApJS 260, 1
``` bibtex
@ARTICLE{2022ApJS..260....1F,
       author = {{Flury}, Sophia R. and {Jaskot}, Anne E. and {Ferguson}, Harry C. and {Worseck}, G{\'a}bor and {Makan}, Kirill and {Chisholm}, John and {Saldana-Lopez}, Alberto and {Schaerer}, Daniel and {McCandliss}, Stephan and {Wang}, Bingjie and {Ford}, N.~M. and {Heckman}, Timothy and {Ji}, Zhiyuan and {Giavalisco}, Mauro and {Amorin}, Ricardo and {Atek}, Hakim and {Blaizot}, Jeremy and {Borthakur}, Sanchayeeta and {Carr}, Cody and {Castellano}, Marco and {Cristiani}, Stefano and {De Barros}, Stephane and {Dickinson}, Mark and {Finkelstein}, Steven L. and {Fleming}, Brian and {Fontanot}, Fabio and {Garel}, Thibault and {Grazian}, Andrea and {Hayes}, Matthew and {Henry}, Alaina and {Mauerhofer}, Valentin and {Micheva}, Genoveva and {Oey}, M.~S. and {Ostlin}, Goran and {Papovich}, Casey and {Pentericci}, Laura and {Ravindranath}, Swara and {Rosdahl}, Joakim and {Rutkowski}, Michael and {Santini}, Paola and {Scarlata}, Claudia and {Teplitz}, Harry and {Thuan}, Trinh and {Trebitsch}, Maxime and {Vanzella}, Eros and {Verhamme}, Anne and {Xu}, Xinfeng},
        title = "{The Low-redshift Lyman Continuum Survey. I. New, Diverse Local Lyman Continuum Emitters}",
      journal = {\apjs},
     keywords = {Reionization, Galactic and extragalactic astronomy, Ultraviolet astronomy, Hubble Space Telescope, 1383, 563, 1736, 761, Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2022,
        month = may,
       volume = {260},
       number = {1},
          eid = {1},
        pages = {1},
          doi = {10.3847/1538-4365/ac5331},
archivePrefix = {arXiv},
       eprint = {2201.11716},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ApJS..260....1F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
## Licensing
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
