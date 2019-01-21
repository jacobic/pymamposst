from collections import OrderedDict
from dataclasses import dataclass
from itertools import chain
from typing import List

from astropy import cosmology

from spiders.objects.misc import WrapperBase


@dataclass
class Data:
    """
    data : str
        data filename is file
    dir : str
        data directory is directory (relative or absolute path)
    DIR_CHAINS : str
        COSMOMC chain directory is directory relative to
        $HOME/NOSAVE/COSMOMC/chains [default: .]
    weightflag : int
        1: Gaussian weights on tracers (with standard deviation taken from
        that of the line-of-sight velocities of the component)
        0: no weights [default: 0]
    distflag : int
        flag for redshift-independent distances [default: 0]:
        0: ignore mu (distance modulus)
        1: weight by Gaussian(mu)
        2: weight by density Gaussian(mu)
    """

    data: str
    dir: str
    chaindir: str = '.'
    weightflag: int = 0
    distflag: int = 0


@dataclass
class Models:
    """
    n : int
        use num tracer components
    components: comp1 [comp2 ... compn]
        use num tracer components with names comp1 ... compn (must be placed
        before other tracer options)
    CofR : file
        [file2 ... filen]
        use files file etc. (with 2 header lines) for R Completeness(R) and
        MAMPOSSt will spline-interpolate [default: none]
    wtofR : file
        [file2 ... filen]
        use files file etc. (with 2 header lines) for R weight(R) and
        MAMPOSSt will spline-interpolate [default: none]
    tracermodel : model
        [model2 ... modeln]
        model(s) model etc. of visible tracer(s); available models are:
        Hernquist: ρ∝r−1(r+a)−3 (Hernquist 1990 model, r−2/a=1/2)
        mHubble: ρ∝(r2+a2)−3/2 (modified Hubble or non-truncated analytical
        King 1962 model, r−2/a=2‾√)
        isothermal: ρ∝(r2+a2)−1 (pseudo-isothermal)
        Jaffe: ρ∝r−2(r+a)−2 (Jaffe 1983)
        NFW: ρ∝r−1(r+a)−2 (Navarro, Frenk & White 1996 model, r−2/a=1)
        Plummer: ρ∝(r2+a2)−5/2 (Plummer 1911 model, r−2/a=2/3‾‾‾√)
        gPlummer: ρ∝rγ(r2+a2)−5/2−γ/2 (generalized Plummer with free inner
        slope γ, r−2/a=(2+γ)/3‾‾‾‾‾‾‾‾√)
        PrugnielSimien: ρ∝x−p(n)exp[−b(n)(rReff)1/n] (Prugniel & Simien 1997
        approximation to deprojected Sersic, where b(n) is from analytical
        approximation of Ciotti & Bertin 1999, while p(n) is given by Lima
        Neto, Gerbal & Marquez 1999)
    anismodel : model
        [model2 ... modeln]
        velocity anisotropy model(s) model1 etc. of visible tracer(s);
        available models are:
        iso: β=0 (isotropic)
        cst: β = constant
        ML: β=12rr+rβ (Mamon & Łokas 2005b)
        OM: β=r2r2+r2β (Osipkov 1979; Merritt 1985)
        gOM: β=β0+(β∞−β0)r2r2+r2β (generalized Osipkov-Merritt)
        Tiret: β=β0+(β∞−β0)rr+rβ (Tiret et al. 2007)
    meanlrtr : meanlr
        [meanlr2 ... meanlrn]
        externally computed log_10 tracer radius (radii) (kpc); tracer radii
        are defined as:
        scale radius a (isothermal or Jaffe)
        effective radius Reff (PrugnielSimien)
        radius of logarithmic density slope −2 (other models, see -tracermodel)
    siglrtr : sigmalr
        [sigmalr2 ... sigmalrn]
        uncertainty of externally computed log_10 tracer radius (radii) (kpc)
    rfidtr : r
        [r2 ... rn]
        tracer fiducial radius (radii) (where tracer mass is given) [kpc,
        default: 1]
    Rminallow : Rmin
        [Rmin2 ... Rminn]
        minimum allowed projected radius (radii) for fits [kpc, 0 for no
        limits, default: 30]
    Rmaxallow : Rmax
        [Rmax2 ... Rmaxn]
        maximum allowed projected radius (radii) for fits [kpc, 0 for no
        limits, default: 0]
    avzmaxallow : vmax
        [vmax2 ... vmaxn]
        maximum allowed | v LOS| for fits (km/s, 0 for no limits, default: 0]
    ev : errv
        [errv2 ... errvn]
        velocity error (in km/s) to be automatically inserted in data file if
        not present [default: not set]

    darkmodel model
       dark matter model; available models are:
       Burkert: ρ∝(r+a)−1(r2+a2)−1 (Burkert 1995 model, r−2/a=(
       1−26/27‾‾‾‾‾√)1/3+(1+26/27‾‾‾‾‾√)1/3)
       Einasto: ρ∝exp[−b(n)r1/n] (Einasto 1965 model, r−2/a=(2n)n)
       Hernquist: ρ∝r−1(r+a)−3 (Hernquist 1990 model, r−2/a=1/2)
       gHernquist: ρ∝rγ(r+a)−4−γ (generalized Hernquist model with free inner
       slope γ, r−2/a=1+γ/2)
       mHubble: ρ∝(r2+a2)−3/2 (modified Hubble or non-truncated analytical
       King 1962 model, r−2/a=2‾√)
       isothermal: ρ∝(r2+a2)−1 (pseudo-isothermal)
       Jaffe: ρ∝r−2(r+a)−2 (Jaffe 1983)
       Kazantzidis: ρ∝rγexp(−r/a) (Kazantzidis et al. 2004 model, r−2/a=2+γ)
       NFW: ρ∝r−1(r+a)−2 (Navarro, Frenk & White 1996 model, r−2/a=1)
       cNFW: ρ∝(r+a)−3 (cored NFW model, r−2/a=2)
       gNFW: ρ∝rγ(r+a)−3−γ (generalized NFW model with free inner slope,
       r−2/a=2+γ)
       Plummer: ρ∝(r2+a2)−5/2 (Plummer 1911 model, r−2/a=2/3‾‾‾√)
       gPlummer: ρ∝rγ(r2+a2)−5/2−γ/2 (generalized Plummer with free inner
       slope γ, r−2/a=(2+γ)/3‾‾‾‾‾‾‾‾√)
    darknormflag dark_norm_flag
       dark matter normalization; available flags are [default: -1]:
       -1: log10(rvir/kpc)
       0: log10(Mvir/M⊙)
       > 0: rfid (kpc) with norm = M(rfid)
    darkscaleflag dark_scale_flag
       dark matter scale flag; available values are [default: 1]:
       1: scale radius
       Reff (PrugnielSimien)
       scale a (isothermal or Jaffe)
       radius of density logarithmic slope −2, r−2 (other models,
       see -darkmodel)
       2: concentration (virial over scale radius)
    darktotflag dark_tot_flag
       dark vs. total flag; available values are [default: 1]:
       1: dark
       2: total
    v3dmodel model
       model for 3D distribution of velocities [default: Gaussian]
    rmax rmax
        max integration radius [kpc, default: 40000]
    """
    components: List
    anismodel: str = None
    CofR = None
    wtofR = None
    tracermodel = None
    meanlrtr = None
    siglrtr = None
    rfidtr = None
    Rminallow = None
    Rmaxallow = None
    avzmaxallow = None
    ev = None
    darkmodel = None
    darknormflag: int = None
    darkscaleflag: int = None
    darktotflag : int = None
    v3dmodel = None
    rmax = None

    def __post_init__(self):
        self.components.insert(0, len(self.components))
        self.n = self.components.copy()
        del self.components


@dataclass
class Variables:
    """
    Variables are in the format
    [start center, min, max, start width, propose width]

     lrtr : min max [min2 max2 ... ]
         log_10 tracer radius (kpc); tracer radii are defined as: scale
         radius a (isothermal or Jaffe) effective radius Reff (
         PrugnielSimien) radius of logarithmic density slope −2 (other
         models, see -tracermodel)
     lMtr : min max [min2 max2 ... ]
         log_10 tracer mass (M⊙) [default: -99 -99 [...]]; tracer masses are
         defined as
         mass within tracer radius (see -lrtr)
     lMtrtot : min max
         log_10 total tracer mass (M⊙) [default: -99 -99]; tracer masses are
         defined
         as mass within tracer radius (see -lrtr)
     ftr tuple min max [min2 max2 ...]
         tracer mass fraction in component (if total tracer mass is set with
     lMtrtot) [default: 0 0 [...]]
     tr2 : tuple
        min max [min2 max2 ...] 2nd tracer parameter (inner slope or Prugniel
        Simien index) [default: 0 0 [...]]
     anisflag : anisotropy_flag
         velocity anisotropy definition flag:
         0: logarithmic: log10(σr/σθ)
         1: standard: β=1−σ2θ/σ2r
         2: symetric: βsym=σ2r−σ2θσ2r+σ2θ
     anis0 : tuple
        min max [min2 max2 ...]
        inner velocity anisotropies (definition in anisotropy_flag) [default:
        0 0 [...]]
     anisinf : tuple
        min max [min2 max2 ...]
         outer velocity anisotropies (definition in anisotropy_flag) [
         default: 0 0 [...]]
     beta0 : tuple
        min max [min2 max2 ...]
         same as -anis0 but force -anisflag 1 [default: 0 0 [...]]
     betainf : tuple
        min max [min2 max2 ...]
         same as -anisinf but force -anisflag 1 [default: 0 0 [...]]
     lranis : tuple
        min max [min2 max2 ...]
         log anisotropy radius (kpc) [default: 0 0 [...]]; can be
         circumvented with
     TALflag :
     norm :  min max
         dark or total (see -darktotflag) mass normalization (see -darknormflag)
     lrdark :
        min max
         log_10 dark or total scale radius (kpc) OR log_10 dark or total
         concentration (depending on dark_scale_flag)
     lc :
        min max
         same as -lrdark min max, but forces -darkscaleflag 2
     darkpar2 :
        min max
         dark matter 2nd parameter (inner slope or Prugniel Simien index) [
         default: 1 1]
     lMBH : min max
         log_10 black hole mass (M⊙) [default: -99 -99]
     lB :
        min max
         log_10 interloper densioty in projected phase space (B parameter in
         Mamon,
         Biviano & Murante 2010, in virial units Nvirr−2vir v −1vir) [
         default: -99 -99]
    """
    lrtr: List = None
    lMtr: List = None
    lMtrtot: List = None
    ftr: List = None
    tr2: List = None
    anisflag: int = None
    anis0: List = None
    anisinf: List = None
    beta0: List = None
    betainf: List = None
    lranis: List = None
    TALflag: int = None
    norm: List = None
    lrdark: List = None
    lc: List = None
    darkpar2: List = None
    lMBH: List = None


@dataclass
class Constraints:
    """
     -TLM
     -TAL
     -split
     -ilopflag
    """

    TLM = None
    TAL = None
    split = None
    ilopflag : int =  None


@dataclass
class Cosmology:
    """
     -z z
    redshift of object [default: 0.0]
    -Omegam Omegam
    cosmological density parameter at z=0, Ω0m [default: 0.3]
    -h h
    Hubble constant at z=0 (units of 100 km/s/Mpc) [default: 0.7]
    -Delta Delta
    value of virial overdensity relative to critical density of Universe at
    observed redshift, 0 for Bryan & Norman 1998 [default: 200]
    -cofM a0 a1
    coefficients for concentratio-nmass relation: log10c=a0+a1log10(M/M⊙) [
    e.g. 2.12 -0.10, default: not set]
    -cofMdef
    assume concentration-mass relation with default coefficients [2.12 -0.10]
    -mu0 mu
    distance modulus of system (if known independently of redshift)
    """

    cosmo: cosmology.core
    z: float = 0.0
    Delta: int = 200

    def __post_init__(self):
        self.Omegam = self.cosmo.Om0
        self.h = (self.cosmo.H0.value / 100)


@dataclass
class MCMC:
    """"
    -samples N
    number N of elements per chain before ending MCMC [default: 104 times the
    number of free ajustable parameters]
    -cv conv
    convergence test on multiple chains: ratio of variance of means to mean
    of variances of chains (e.g. 1, 0.1, or 0.01; default: not set)
    -p num
    propose scale [default: 2.4]
    -burn_in N
    wait for N “burn in” elements per chain before speeding up search using
    covariance matrix calculations [default: 1000]
    -feedback f
    CosmoMC verbosity [default: 0]
    -seed integer
    random seed: 0 for seed based on clock [default: 1234]
    """

    samples: int = None
    cv: float = None
    p: float = 2.4
    burn_in: int = 1000
    feedback: int = 0
    seed: int = 1234


@dataclass
class Extra:
    """

    Parameters
    ----------
    MAMPOSStParams

    Returns
    -------

    """

    prefix: str
    v = None
    debug: int = 0


class Builder(WrapperBase):
    """

    Build MAMPOSSt initialistion file

    buildini_mamposst - data_dir. - data_file test.dat - DIR_CHAINS TESTS - n 2
    Passive SF - lrtr 2.5 3.4 3 3.9 test
    Returns
    -------

    """

    def __init__(self, data, models, extra, variables=None, constraints=None,
                 cosmology=None, mcmc=None):
        self.data = data
        self.cosmology = cosmology
        self.models = models
        self.variables = variables
        self.constraints = constraints
        self.mcmc = mcmc
        self.extra = extra
        arg_groups = [data, cosmology, constraints, models, variables, mcmc,
                      extra]

        # Must preserve order during update of dictionary with OrderedDict.
        self._args = OrderedDict({})
        for group in arg_groups:
            if group is not None:
                group_dict = OrderedDict(group.__dict__)
                if isinstance(group, Models):
                    # Enter -n ncomp comp1 .. before tracer quantities or models
                    if 'n' in group_dict.keys():
                        group_dict.move_to_end('n', last=False)
                elif isinstance(group, Cosmology):
                    # Remove cosmology.core object to leave cmdline params.
                    group_dict.pop('cosmo')
                self._args.update(group_dict)

        self.args = {'-{0}'.format(k): v for k, v in self._args.items()}

        def build_args(arg_dict):
            for k, v in arg_dict.items():
                if v is not None:
                    # Value is not None.
                    if isinstance(v, list):
                        # Value is actually a list of values.
                        arg_list = [k.split(), [str(value) for value in v]]
                        arg = list(chain.from_iterable(arg_list))
                    else:
                        # Value is a int or str
                        arg = [k, str(v)]
                    yield arg

        cmd_args = list(chain.from_iterable(list(build_args(self.args))))
        cmd = ['buildini_mamposst'] + cmd_args
        # -prefix option does not actually exist, it is just for convenience.
        cmd.remove('-prefix')
        super().__init__(cmd=cmd, logname='_build.log')

    def build(self):
        """
        buildini_mamposst -data_dir. -data_file test.dat -DIR_CHAINS TESTS -n 2
        Passive SF -lrtr 2.5 3.4 3 3.9 test

        Returns
        -------

        """
        self.run_cmd()
