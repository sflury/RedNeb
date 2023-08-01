'''
dependencies:
    numpy 1.17.3
    pyneb 1.1.10
    scipy 1.4.1
    pandas 1.2.1

note:
    pandas dependency for dered_iter only

author:
    Sophia Flury 2022.02.15
'''
import numpy as np
import pandas as pd
import pyneb as pn
from scipy.optimize import curve_fit
# set up and instantiate PyNeb recombination object
global H1
H1 = pn.RecAtom('H',1)


'''
Function to artificially redden intrinsic Balmer decrements and introduce
stellar absorption effects. Used by calc_ebmv_from_balmer to fit intrinsic
decrements obtained from temperatures and densities to observed decrements.

Model is
    f_obs/f_obs,H-beta = f_intr/f_intr,Hbeta * (EW_obs + EW_abs)/(EW_obs) * \
                          (EW_obs,Hbeta)/(EW_obs,Hbeta + EW_abs) * \
                          10**( -0.4*E[B-V]*(k_H-k_Hbeta) )

'''
def balmer_red(waves,f_intr,ew,ebmv,star_abs,ext_law='ccm89',j=int(0)):
    # add stellar absorption by scaling decrement
    #ystel = x*EWe/(EWe+ew_a)*(EWeb+ew_a)/EWeb
    f_star_corr = f_intr*ew/(ew+star_abs)*(ew[j]+star_abs)/ew[j]
    # add reddening
    f_red_corr = f_star_corr * 10**(-0.4*ebmv*get_ext_coeff(waves,ext_law))
    return f_red_corr/f_red_corr[j]

'''
Nebular functions taken from nebular_classes.py to remove object class component
'''
### get OII temp from OIII temp using Perez-Montero and Diaz 2003 MNRAS 346, 105
def o2temp(Te,ne):
    return 1.e4/(0.693/(Te/1e4)+0.281)
### get OII temp from OIII temp using Dors et al 2020 MNRAS XXX, XXX
def o2temp_dors20(Te,ne):
    return 1.e4*np.polyval([0.17,-1.07,2.07,-0.33],Te/1e4)

### iterative solving of equilibrium equations using PyNeb routines
def calc_tempdens(t_ratio,d_ratio,t_atom='O3',d_atom='S2',n_iter=int(15),tol=2**-5):
    # dict of PyNeb Atom objects for TemDen diagnostics
    atoms = {}
    for atom in [t_atom,d_atom]:
        if atom not in atoms.keys():
            atoms[atom] = pn.Atom(atom[:-1],float(atom[-1]))
    # check if array or list of ratios
    if hasattr(t_ratio,'__len__'):
        # run converging solver for each set of ratios
        return np.array(list(map(lambda tr,dr: \
                __tempdens_converge(tr,dr,atoms,n_iter=n_iter,tol=tol),\
                t_ratio,d_ratio))).T
    # if not array, just run once
    else:
        return __tempdens_converge(t_ratio,d_ratio,atoms, \
                    n_iter=n_iter,tol=tol)

def __tempdens_converge(t_ratio,d_ratio,atoms,n_iter=int(15),tol=1.):
    # t_atom from atoms dict
    t_atom = list(atoms.keys())[0]
    # d_atom from atoms dict
    if len(atoms.keys()) == 2:
        d_atom = list(atoms.keys())[1]
    else:
        d_atom = list(atoms.keys())[0]
    # dicts of supported levels
    tlvls = {'S2':'(I(5,1)+I(5,2))/(I(3,1)+I(2,1))', \
                    'O2':'(I(3,1)+I(2,1))/(I(4,2)+I(4,3)+I(5,2)+I(5,3))',\
                    'N2':'I(5,4)/(I(4,2)+I(4,3))',\
                    'O3':'I(5,4)/(I(4,2)+I(4,3))'}
    dlvls = {'S2':'L(6731)/L(6716)','O2':'L(3726)/L(3729)'}
    # check on density ratios -- if too low, set to low-density limit
    if d_ratio < atoms[d_atom].getLowDensRatio(to_eval=dlvls[d_atom]):
        d_ratio = 1.1*atoms[d_atom].getLowDensRatio(to_eval=dlvls[d_atom])
    # initial conditions
    Ti = atoms[t_atom].getTemDen(t_ratio,den=100,to_eval=tlvls[t_atom])
    ni = atoms[d_atom].getTemDen(d_ratio, \
                    tem=o2temp(Ti,100),to_eval=dlvls[d_atom])
    # iteratively solve for Te and ne until converging
    # or max iteration count is reached
    n = int(1)
    while n <= n_iter :
        # estimate temperature with density guess
        Te = atoms[t_atom].getTemDen(t_ratio,den=ni,to_eval=tlvls[t_atom])
        # estimate density with temperature guess
        ne = atoms[d_atom].getTemDen(d_ratio, \
                    tem=o2temp(Ti,ni),to_eval=dlvls[d_atom])
        # NaN check on density
        if np.isnan(ne) :
            ne = 100.
        # NaN check on temperature
        if np.isnan(Te) :
            Te = 1e4
        # if converged on solution, return
        if abs(Te-Ti) < tol and abs(ne-ni) < tol :
            return Te,ne
        # otherwise, update and continue
        else:
            Ti = np.copy(Te)
            ni = np.copy(ne)
        n += 1
    return Te,ne

'''
Purpose
    Obtain extinction coefficients for desired extinction law
    at given wavelengths.

Arguments
    :wave (*float* or *np.ndarray*): wavelength(s) in Ang at which to evaluate
            the extinction law

Keyword Arguments
    :ext_law (*str*): Extinction law to evaluate. Currently
            supports \'ccm89\', \'calzetti\',\'smc\', and
            \'fitz99\'. Default is \'ccm89\'.

Returns
    Extinction coefficient(s) A_lambda / E(B-V) evaluated at the given
    wavelength(s) for the desired extinction law
'''
def get_ext_coeff(wave,ext_law='ccm89'):
    w_ext = [3655,3889,3970,4102,4341,4861,6563]
    alamb_ebmv_ext = \
            {'ccm89':[4.73821,4.55620,4.48243,4.35016,4.08539,3.51976,2.45495],\
            'calzetti':[6.78098,6.47902,6.37941,6.22219,5.95279,5.42747,4.19790],
            'smc':[4.51798,4.25315,4.15807,4.00515,3.73772,3.19774,2.11598],\
            'fitz99':[4.68657,4.45715,4.38187,4.26244,4.05112,3.55609,2.30253]}
    return np.interp(wave,w_ext,alamb_ebmv_ext[ext_law])

'''
Purpose
    Apply reddening correction to given fluxes at given wavelengths
    for a value of E(B-V) and desired extinction law.

Arguments
    :wave (*float* or *np.ndarray*): wavelength(s) in Ang at which to evaluate
            the extinction law
    :flux (*float* or *np.ndarray*): fluxes corresponding to wave which to
            correct for reddening
    :eBmV (*float*): E(B-V) value to scale extinction coefficients for reddening
            correction

Keyword Arguments
    :ext_law (*str*): Extinction law to evaluate. Currently
            supports \'ccm89\', \'calzetti\',\'smc\', and
            \'fitz99\'. Default is \'ccm89\'.

Returns
    Reddening-corrected fluxes for the assumed extinction law scaled by given
    E(B-V).
'''
def red_corr(wave,flux,eBmV,ext_law='ccm89'):
    return flux * 10**(0.4*eBmV*get_ext_coeff(wave,ext_law))

'''
Name
    calc_ebmv_from_balmer

Purpose
    Calculate E(B-V) for a desired extinction using Balmer line fluxes
    alpha, beta, gamma, and delta.  Computes recombination coefficients
    and extinction coefficients to populate a NxN matrix containing
    E(B-V) values from all N available Balmer fluxes. Rows and columns
    organized as below by the pair of lines used to compute E(B-V) where
    the ith row is the numerator and the jth column is the denominator
    (although this does not matter since the flux ratio cancels the sign).
    Diagonal elements set to nan since a single line yields no
    selective extinction. Example matrix if Ha, Hb, Hy, and Hd fluxes
    provided to the script:

                  Ha  Hb  Hy  Hd
               Ha --
               Hb     --
               Hy         --
               Hd             --

    Returns the upper diagonal elements of this matrix, ordered row by row.

Arguments
    :wave (*np.ndarray*): 1xN array of the Balmer line wavelengths
    :flux (*np.ndarray*): 1xN array of the Balmer line fluxes

Keyword Arguments
    :tem (*float*): Balmer electron temperature in K. Default is 10^4 K.
    :den (*float*): Balmer electron density in cm^-3. Default is 10^2 cm^-3.
    :ext_law (*str*): String containing the name of the extinction law.
            Currently supports \'ccm89\', \'calzetti\',\'smc\', and
            \'fitz99\'. Default is \'ccm89\'.

Returns
    :eBmV (*np.ndarray*): a 1xC array of E(B-V) values derived from the
            observed Balmer flux ratios using Case B recombination
            coefficients and an assumed extinction law. When flux is not
            reported or computed E(B-V) is negative, E(B-V) is set to zero.
            C is given by the number of combinations of flux ratios without
            duplicates such that C=1/2 N!/(N-2)! for n emission line fluxes.
'''
def calc_ebmv_from_balmer(wave,flux,flux_err=None,temp=1e4,dens=1e2,ext_law='ccm89'):
    # number of lines
    n = len(flux)
    # starting energy levels
    B = 3645.0682 # Balmer's constant
    Narr = [int(round(2*np.sqrt(np.max(w)/(np.max(w)-B)))) for w in wave]
    # extinction coefficients by interpolating over extinction law
    k_lamb = get_ext_coeff(wave,ext_law=ext_law)
    # recombination coefficients for alpha through delta
    alpha = np.array([H1.getEmissivity(tem=temp,den=dens,wave=w) for w in wave])
    # array of E(B-V) values
    eBmV = np.zeros((n,n))
    eBmV_err = np.zeros((n,n))
    # for all lines
    for i in range(n):
        # for all possible lines (0 to i-1 are duplicates,
        # i==j has no selective extinction)
        for j in range(i+1,n):
            # get E[B-V] from decrements and extinction coefficients
            # following Fi/Fj = F0i/F0j * 10^(-0.4E[B-V](k_lamb,i-k_lamb,j))
            eBmV[i,j] = -2.5/(k_lamb[i]-k_lamb[j]) * \
                        np.log10((flux[i]/flux[j])/(alpha[i]/alpha[j]) )
            # get error in E[B-V]
            if not np.any(flux_err==None):
                eBmV_err[i,j] = 2.5/abs(k_lamb[i]-k_lamb[j])/np.log(10) * \
                        np.sqrt( (flux_err[i]/(flux[i]))**2 + \
                        (flux_err[j]/(flux[j]))**2)
    # return upper diagonal elements
    inds = np.triu_indices(n,k=1)
    return eBmV[inds],eBmV_err[inds]

'''
Name
    fit_ebmv_starabs

Purpose
    Calculate E(B-V) and stellar absorption correction for a desired extinction
    law using provided Balmer line fluxes.  Computes recombination coefficients
    using PyNeb and extinction coefficients using the assumed extinction law.
    Uses scipy.optimize.curve_fit implementation of Levenberg-Marquardt to
    fit intrinsic emissivities to observed fluxes with E(B-V) and stellar
    absorption lines as free parameters.

Arguments
    :flux (*np.ndarray*): 1xn array of the Balmer line fluxes ordered by
            decreasing wavelength as floats, starting with H alpha
            (3->2, w=6563 Ang), proceeding to H beta (4->2, w=4861 Ang),
            H gamma (5->2, w=4340 Ang), etc. Substitute zeros for
            undetected lines. Supports up to H eta (9->2, w=3835 Ang).
    :ew (*np.ndarray*): 1xn array of the Balmer line equivalent widths ordered
            in the same manner as `flux`.

Keyword Arguments
    :ebv_guess (*float*): Best guess for E(B-V). Default is None.
    :tem (*float*): Balmer electron temperature in K. Default is 10^4 K.
    :den (*float*): Balmer electron density in cm^-3. Default is 10^2 cm^-3.
    :ext_law (*str*): String containing the name of the extinction law.
            Currently supports \'ccm89\', \'calzetti\',\'smc\', and
            \'fitz99\'. Default is \'ccm89\'.

Returns
    :eBmV (*np.ndarray*): best-fit E(B-V), accounting for uncertainty in Balmer
                        lines and equivalent widths
    :star_abs (*np.ndarray*): best-fit stellar absorption in Ang, accounting for
                        uncertainty in Balmer lines and equivalent widths
    :eBmV_err (*np.ndarray*): uncertainty in best-fit E(B-V)
    :star_abs_err (*np.ndarray*): uncertainty in best-fit stellar absorption

'''
def fit_ebmv_starabs(wave,flux,ew,flux_err=None,ebmv_guess=None,temp=1e4,dens=1e2,ext_law='ccm89'):
    #clp = np.where(np.isfinite(flux))
    #wave = wave[clp]
    #flux = flux[clp]
    #ew   = ew[clp]
    #if not np.any(flux_err == None): flux_err = flux_err[clp]
    # starting energy level
    B = 3645.0682 # Balmer's constant in Ang
    N = [int(round(2*np.sqrt(w/(w-B)))) for w in wave]
    # recombination coefficients
    alpha = np.array([H1.getEmissivity(tem=temp,den=dens,wave=w) for w in wave])
    if ebmv_guess == None:
        # E(B-V) from emissivities
        ebmv_arr,ebmv_errs = calc_ebmv_from_balmer(wave,flux,\
                            flux_err=flux_err,ext_law=ext_law)
        weight = ebmv_errs**-2
        ebmv_guess = np.nansum(ebmv_arr*weight)/np.nansum(weight)
        ebmv_guess_err = np.sqrt(1./np.nansum(weight))
    # locate H-beta
    j = np.argmin(abs(wave-4861))
    # if no h-beta, just use zeroth line to normalize
    if abs(wave[j]-4861) > 2. :
        j = int(0)
    # boundary conditions and initial conditions based on initial guess
    if ebmv_guess > 0:
        bounds = [[0.0,0.0],[2.,min([50,ew[j]])]]
        p0 = [ebmv_guess,0.0]
    # if initial guess is negative (on average, decrements forbidden)
    else:
        bounds = [[0.0,0.0],[2.,min([50,ew[j]])]]
        p0 = [0.05,0.0]
    # hard-wire equivalent widths
    func = lambda f,ebmv,star_abs: balmer_red(wave,f,ew,ebmv,star_abs,\
                        ext_law=ext_law,j=j)
    # LM fit to get E(B-V), accounting for errors
    par_opt,covar = curve_fit(func,alpha/alpha[j],flux/flux[j], \
                        sigma=flux_err,p0=p0,bounds=bounds,maxfev=10000)
    sigma = np.sqrt(np.diag(covar))
    return par_opt[0],par_opt[1],sigma[0],sigma[1]

'''
Name
    fit_starabs

Purpose
    Calculate stellar absorption correction for a desired extinction law and
    known/assumed E(B-V) using provided Balmer line fluxes and equivalent
    widths.  Computes recombination coefficients using PyNeb and extinction
    coefficients using the assumed extinction law.Uses scipy.optimize.curve_fit
    implementation of Levenberg-Marquardt to fit intrinsic emissivities to
    observed fluxes with stellar absorption lines as the free parameter.

Arguments
    :flux (*np.ndarray*): 1xn array of the Balmer line fluxes ordered by
            decreasing wavelength as floats, starting with H alpha
            (3->2, w=6563 Ang), proceeding to H beta (4->2, w=4861 Ang),
            H gamma (5->2, w=4340 Ang), etc. Substitute zeros for
            undetected lines. Supports up to H eta (9->2, w=3835 Ang).
    :ew (*np.ndarray*): 1xn array of the Balmer line equivalent widths ordered
            in the same manner as `flux`.

Keyword Arguments
    :ebmv (*float*): Assumed E(B-V). Default is 0.0.
    :tem (*float*): Balmer electron temperature in K. Default is 10^4 K.
    :den (*float*): Balmer electron density in cm^-3. Default is 10^2 cm^-3.
    :ext_law (*str*): String containing the name of the extinction law.
            Currently supports \'ccm89\', \'calzetti\',\'smc\', and
            \'fitz99\'. Default is \'ccm89\'.

Returns
    :star_abs (*np.ndarray*): best-fit stellar absorption in Ang, accounting for
                        uncertainty in Balmer lines and equivalent widths
    :star_abs_err (*np.ndarray*): uncertainty in best-fit stellar absorption

'''
def fit_starabs(wave,flux,ew,ebmv=0.0,flux_err=None,temp=1e4,dens=1e2,ext_law='ccm89'):
    # recombination coefficients
    alpha = np.array([H1.getEmissivity(tem=temp,den=dens,wave=w) for w in wave])
    # locate H-beta
    hb = np.argmin(abs(wave-4861))
    # boundary conditions and initial conditions based on initial guess
    bounds = [[0.0],[20.]]
    p0 = [0.0]
    # hard-wire equivalent widths
    func = lambda f,star_abs: balmer_red(f,ew,ebmv,star_abs,ext_law=ext_law)
    # LM fit to get E(B-V), accounting for errors
    par_opt,covar = curve_fit(func,alpha/alpha[hb],flux/flux[hb], \
                        sigma=flux_err,p0=p0,bounds=bounds,maxfev=10000)
    sigma = np.sqrt(np.diag(covar))
    return par_opt[0],sigma[0]

'''
Name
    dered_iter

Purpose
    Compute E(B-V) from Balmer lines using the temperature and density computed
    from the emission lines to determine the intrinsic Balmer decrement. Fluxes
    are corrected for inferred, error-weighted E(B-V) value, temp and dens are
    recomputed from corrected lines, and process is iteratively repeated until
    converging on a consistent solution for E(B-V).

Arguments
    :flux0 (*pd.DataFrame*): single-row pandas dataframe with emission line
                        fluxes anderrors in columns labeled in PyNeb fashion:
                        ElemSpec_waveA and ElemSpec_waveAe. Recombination lines
                        require additional r after the species number. E.g., for
                        H-alpha, H1r_6563A for flux and H1r_6563Ae for error and
                        for [O III]5007, O3_5007A for flux and O3_5007Ae for
                        error. Fluxes are assumed to have been corrected for
                        MilkyWay extinction prior to passing to dered_iter.
    :eqwd (*pd.DataFrame*): single-row pandas dataframe with emission line
                        fluxes and errors in columns labeled as `flux0`.

Keyword Arguments
    :ext_law (*str*): string indicating the extinction law to use. Currently
                        supports \'ccm89\', \'calzetti\',\'smc\', and
                        \'fitz99\'. Default is \'ccm89\'.
    :abs_corr (*bool*): boolean indicating whether to include stellar absorption
                        correction in the E(B-V) estimate. If `True`, will use
                        scipy.optimize.curve_fit to fit the intrinsic ratios to
                        the observed ratios of Balmer fluxes assuming

                            f_obs/f_obs,Hb = f_0/f_0,Hb*10^(-0.4*E[B-V]*(k-k_Hb))
                                            *(ew+ew_abs)/ew*ew_Hb/(ew_Hb+ew_abs)

                        If `False`, will calculate the E(B-V) for all possible
                        Balmer flux ratio permutations assuming

                            f_obs1/f_obs2 = f_0,1/f_0,2*10^(-0.4*E[B-V]*(k_1-k_2))

                        and take the error-weighted average to be the
                        characteristic extinction.
                        Default is `False`.
    :verbose (*bool*): boolean indicating whether to print used Balmer lines,
                        observed fluxes (normalized to H-beta), corrected fluxes
                        (normalized to H-beta), and characteristic E(B-V) and
                        stellar absorption correction.

Returns
    :ebmv (*float*): variance-weighted E(B-V) estimated from Balmer lines
    :ebmv (*float*): error in E(B-V). If `abs_corr=True`, this is computed from
                        the covariance matrix of the best-fit parameters. If
                        `abs_corr=False` (default), this is computed from the
                        variance weights used to compute the characteristic
                        E(B-V).
    :star_abs (*float*): if `abs_corr=True`, the error-weighted stellar
                        absorption correction; if `abs_corr=False`, zero.
    :star_abs (*float*): if `abs_corr=True`, the error in the stellar
                        absorption correction computed from the covariance
                        matrix of best-fit parameters; if `abs_corr=False`,
                        zero.
'''
def dered_iter(flux0,eqwd,ext_law='ccm89',abs_corr=False,verbose=True):
    # flux dict from pandas dataframe object
    lines = np.array([key for key in flux0.keys()[1:] if not key.endswith('e')])
    # wavelengths from column headers
    wave = np.array([float(l.split('_')[1][:-1]) for l in lines])
    # find measured balmer fluxes
    balmer = np.array([ i for i,l in enumerate(lines) if l.startswith('H1r_') \
                        if flux0[l].values[0]/flux0[l+'e'].values[0]>3 ])[::-1]
    if len(balmer) < 2 :
        if verbose:
            print('Insufficient number of Balmer lines to derive extinction.')
        return 0.,0.,0.,0.
    # get uncorrected balmer fluxes
    f_balmer = np.array([flux0[lines[j]].values[0] for j in balmer])
    f_balmer_err = np.array([flux0[lines[j]+'e'].values[0] for j in balmer])
    ew = np.array([eqwd[lines[j]].values[0] for j in balmer])
    # find H-alpha and H-beta
    ind_ha = np.argmin(abs(wave[balmer]-6563))
    ind_hb = np.argmin(abs(wave[balmer]-4861))
    # check for bad H-alpha/H-beta and remove unphysical ratios
    hahb = f_balmer[ind_ha]/f_balmer[ind_hb]
    hahbe = np.sqrt((f_balmer_err[ind_ha]/f_balmer[ind_ha])**2.+\
                        (f_balmer_err[ind_hb]/f_balmer[ind_hb])**2.)
    if hahb + hahbe < 2.74 and len(f_balmer)>2:
        balmer = balmer[1:]
        f_balmer = f_balmer[1:]
        f_balmer_err = f_balmer_err[1:]
        ew = ew[1:]
        ind_hb -= 1
        if verbose:
            print('Unphysical j6563/j4861. H-alpha not used.')
    # detected lines
    if verbose:
        print('Detected Balmer lines')
        ltxt = (len(f_balmer)*'\tj{:d}/j4861')
        print(ltxt.format(*[int(wave[b]) for b in balmer]))
        print('obsv:'+(len(f_balmer)*'\t{:f}').format(*f_balmer/f_balmer[ind_hb]))
    # new flux dict
    flux = {}
    for l in lines:
        flux[l] = flux0[l].values
    # check if not yet converged
    check = True
    # start with no reddening
    ebmv0 = 0.
    # counter
    n_iter = int(0)
    # iterate until convergence
    while check and n_iter < 20 :
        # compute temperature and density
        t_ratio = flux['O3_4363A']/(flux['O3_4959A']+flux['O3_5007A'])
        d_ratio = flux['S2_6731A']/flux['S2_6716A']
        te,ne = calc_tempdens(t_ratio[0],d_ratio[0])
        # get E(B-V) from Balmer lines
        # with correcting for stellar absorption
        if abs_corr:
            # Levenberg-Marquardt fit to Balmer fluxes and equivalent widths
            ebmv1,star_abs,ebmv_err,star_abs_err = \
                fit_ebmv_starabs(wave[balmer],f_balmer,ew,\
                flux_err=f_balmer_err,ext_law=ext_law,
                temp=te,dens=ne)
        # without correcting for stellar absorption
        else:
            # get E(B-V) from all possible Balmer flux ratios
            ebmv_arr,ebmv_errs = calc_ebmv_from_balmer(wave[balmer],f_balmer,\
                            flux_err=f_balmer_err,ext_law=ext_law,
                            temp=te,dens=ne)
            # error-weighted average
            weight = ebmv_errs**-2
            ebmv1 = np.nansum(ebmv_arr*weight)/np.nansum(weight)
            ebmv_err = np.sqrt(1./np.nansum(weight))
            # dummy star absorption variables
            star_abs = 0.0
            star_abs_err = 0.0
        # compare old and new E(B-V)
        check = abs(ebmv1-ebmv0) > 2**-16
        # update comparison extinction
        ebmv0 = np.copy(ebmv1)
        # update counter
        n_iter += 1
        # if not yet updated, correct fluxes using best-fit extinction
        for i,l in enumerate(lines):
            # correct fluxes for extinction
            flux[l] = red_corr(wave[i],flux0[l].values,ebmv1,ext_law=ext_law)
            # correct Balmer fluxes for stellar absorption
            if l.startswith('H1r_'):
                flux[l] *= (eqwd[l].values[0]+star_abs)/eqwd[l].values[0]
    # outputs
    f_balmer = np.array([flux[lines[j]][0] for j in balmer])
    # results
    if verbose:
        ctxt = 'corr:'+(len(f_balmer)*'\t{:f}')
        print(ctxt.format(*f_balmer/f_balmer[ind_hb]))
        rtxt = '\tN_ITER: {:d}\tE(B-V): {:2.6f}\tStellar Abs EW: {:2.6f}'
        print(rtxt.format(n_iter,ebmv0,star_abs))
    # return converged value
    return ebmv1,ebmv_err,star_abs,star_abs_err
