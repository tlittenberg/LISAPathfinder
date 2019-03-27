import numpy as np
import matplotlib.pyplot as plt
from populationClass import population as pop
from impactClass import impactClass
from impactClass import impactClassList
import pathlib
import os 
import re
import time
import pandas as pd
import sys


p = pathlib.PurePath(os.getcwd())
BASE_DIR = str(p.parent)
dataDir = '/data'

# Number of GOOD impacts, should save this number somewhere
N_DETECTED = 44

def getGRSSegments(filenames, grs = 1):
        # Get list of segment times
        regex = r'(\d*)_grs%i.*'%(grs)
        segments = []
        for f in filenames:
                if re.match(regex, f):
                        segments.append(int(re.findall(regex, f)[0]))
        return segments


def redefine(dataDir = '/data', saveDir = '/FLUX_IMPACTS_ALL',usePtot=True):
        """
        Creates .npy files used in the analysis

        """
        dataPath = pathlib.Path(BASE_DIR + dataDir)
        #usePtot = True
        if not usePtot and not '_INTP' in saveDir:
                sys.exit('Expect dir name FLUX_IMPACTS_ALL_INTP with usePtot=False')
        
        # Set up directory structure
        impact_dir = pathlib.Path(str(dataPath) + '/ALL_IMPACTS')
        save_dir = pathlib.Path(str(dataPath) + saveDir)
        modelDir = pathlib.Path(str(dataPath) + '/models')

        # Read in populations
        print("Initializing Populations")
        populations = [ pop(modelDir, 'JFC', usePtot),
                                        pop(modelDir, 'HTC', usePtot),
                                        pop(modelDir, 'AST', usePtot),
                                        pop(modelDir, 'OCC', usePtot),
                                        pop(modelDir, 'Uniform', usePtot)]

        pickles = list(impact_dir.glob('*_grs1.pickle'))

        print("Reading through pickle files")
        grs = 1

        # Set up likelihood file
        # Write likihood file header
        write_filename = str(save_dir) + '/rename.log'
        write_file = open(write_filename, 'w+')
        write_file.write('Making numpy files')
        write_file.write('\n')

        i = -1
        for p in pickles:
                        i += 1
                        segment = str(p.stem[0:10])
                        print(segment, i * 100 / len(pickles), '%')
                        write_file.write(segment + '\n')
                        file_numpy = str(save_dir) + '/' + segment + '_grs%i.npy'%(grs)

                        if os.path.isfile(file_numpy):
                                print("Exists")
                                continue
                        chainFile = str(impact_dir) + '/' + str(segment) +'_grs%i'%grs + '.pickle'

                        # Initialize Impact
                        impact = impactClass(chainFile)
                        impact.segment = int(segment)
                        impact = impact.SCtoSun()
                        impact = impact.findSkyAngles()

                        if impact.lon_sun is None:   #JGB: This indicates no impacts at all in the chain.  There are 420 such examples, verified in a few examples.
                                print( "NONE")
                                lons = [None]
                                lats = [None]
                                flux_JFC = [None]
                                flux_HTC = [None]
                                flux_AST = [None]
                                flux_OCC = [None]
                                flux_Uniform = [None]
                        else:
                                try:
                                        lons = np.asarray(impact.lon_sun)
                                        lats = np.asarray(impact.lat_sun)
                                        Ptots = np.asarray(impact.Ptot)
                                        if usePtot:
                                                flux_JFC = populations[0].getFlux(lons, lats, Ptots, norm = False)
                                                flux_HTC = populations[1].getFlux(lons, lats, Ptots,  norm = False)
                                                flux_AST = populations[2].getFlux(lons, lats, Ptots, norm = False)
                                                flux_OCC = populations[3].getFlux(lons, lats, Ptots,  norm = False)
                                                # This number is meaningless with no norm but we dont have it for simplification reading in 
                                                flux_Uniform = populations[4].getFlux(lons, lats, Ptots,  norm = False)
                                        else:
                                                flux_JFC = populations[0].getFlux(lons, lats, norm = False)
                                                flux_HTC = populations[1].getFlux(lons, lats, norm = False)
                                                flux_AST = populations[2].getFlux(lons, lats, norm = False)
                                                flux_OCC = populations[3].getFlux(lons, lats, norm = False)
                                                # This number is meaningless with no norm but we dont have it for simplification reading in 
                                                flux_Uniform = populations[4].getFlux(lons, lats, norm = False)

                                #except:
                                except Exception as e:
                                        print("Something went wrong: ",str(e))
                                        print("Trying again")
                                        fluxes = np.zeros((len(lons), len(populations)))
                                        for i, population in enumerate(populations):
                                                for j in range(len(lons)):
                                                        fluxes[j, i] = populations[i].getFlux(lons[j], lats[j], Ptots[j], norm = False)
                                        flux_JFC = fluxes[:, 0]
                                        flux_HTC = fluxes[:, 1]
                                        flux_AST = fluxes[:, 2]
                                        flux_OCC = fluxes[:, 3]
                                        flux_Uniform = fluxes[:, 4]

                        
                        print(impact.segment,lons[0:3],lats[0:3],flux_JFC[0:3],flux_HTC[0:3],flux_Uniform[0:3])
                        N = impact.N * np.ones(len(lons))
                        segment = impact.segment * np.ones(len(lons))

                        data_array = [impact.lon_sun, impact.lat_sun, impact.Ptot, 
                                                  flux_JFC, flux_HTC, flux_AST, flux_OCC, flux_Uniform, N, segment]

                        np.save(file_numpy, data_array)
                        write_file.write('saved' + '\n')
                        write_file.write('\n')

class simpleImpact:
        def __init__(self, segment, direc = '/data/FLUX_IMPACTS_ALL',usePtot=True):
                segment_direc = BASE_DIR + direc
                if not usePtot:
                        segment_direc+='_INTP'
                data = np.load(segment_direc + '/' + str(segment) + '_grs1.npy')
                #print("reading",segment)
                #print(data.shape)
                #print(np.array(data))
                try:
                        self.lon = data[0, :]
                        self.lat = data[1, :]
                        self.Ptot = data[2, :]
                        # JFC, HTC, AST, OCC, Uniform
                        self.flux = data[3:8, :]

                        self.N = int(data[8, 0])
                        self.Nimpact=len(self.lon)
                        self.segment = int(data[9, 0])

                except:  #JGB:  Indicates no impacts in chain
                        #print( "None case:",segment)
                        self.lon = [None]
                        self.lat = [None]
                        self.Ptot = [None]
                        self.flux = [None]
                        self.N = int(data[8])
                        self.Nimpact=0
                        self.segment = int(data[9])

                #print('lon:',self.lon)
                #print('lat:',self.lat)
                #print('Ptot:',self.Ptot)
                #print('flux:',self.flux)
                #print('N/Ni:',self.N,self.Nimpact)
                #print('seg:',self.segment)
                #print('flux[:,0:3]=',self.flux[:,0:3])
                #print('flux[:,1000:1003]=',self.flux[:,1000:1003])
                return

def getVetoList():
        df_veto = pd.read_csv(BASE_DIR + dataDir + '/impact_list.txt', header = 'infer',
                        delim_whitespace = True)
        return df_veto


def getisGlitch(segment, df_veto):
        try:
                index = df_veto.index[df_veto['segment'] == int(segment)][0]
        except IndexError:
                # if not in the list, its not an Impact
                return False
        else:
                return df_veto['isGlitch'].values[index]

def getisValidSearch(segment):
        search_times = pd.read_csv(BASE_DIR + dataDir + '/segment_list.txt',
                        header = None, names = ['segment'], delim_whitespace = True)

        # Jakes list does not cover all time, keep segments only within list
        if segment > max(search_times['segment'].values):
                return True
        elif int(segment) in search_times['segment'].values:
                return True
        else:
                return False

def isValid(segment, df_veto):
        #       Both checks if the time searched is valid as well as 
        #       if the time is a glitch

        if getisValidSearch(segment) and not getisGlitch(segment, df_veto):
                return True

        # If it's not in the segment list then its not a good segment
        return False

def getIndex(populations):
        # This is the order of the populations in the 
        # .npy files
        index = []
        for p in populations:
                if p.pop_type == 'JFC':
                        index.append(0)
                elif p.pop_type == 'HTC':
                        index.append(1)
                elif p.pop_type == 'AST':
                        index.append(2)
                elif p.pop_type == 'OCC':
                        index.append(3)
                elif p.pop_type == 'Uniform':
                        index.append(4)
        return index


def calcPQuantiles(segments,Nquant=1000):
        print("counting all momenta")
        #First get len of smallest segment
        Nmin=1e100;
        for segment in segments:
                Nmin=min([Nmin,segment.N])
        #Next get total # of momenta after scaling segments to same size
        NP=0
        for segment in segments:
                if segment.Ptot[0] is not None: 
                        trimNp=int(len(segment.Ptot)*Nmin/float(segment.N)+0.5)
                        NP+=trimNp
        allN=len(segments)*Nmin;
        print('allN,Nmin,NP',allN,Nmin,NP)
        #Next make a vector of all momenta (it is way faster to do this in two steps, to avoid repeatedly resizing the array)
        allP=np.zeros(NP)
        i=0;
        for segment in segments:
                if segment.Ptot[0] is not None: 
                        trimNp=int(len(segment.Ptot)*Nmin/float(segment.N)+0.5)
                        trimP=np.random.choice(segment.Ptot,trimNp)
                        allP[i:i+trimNp]=trimP
                        i+=trimNp                        
        Nunit=Nmin
        xs=np.array([x for x in np.arange(0,NP,NP/float(Nquant))])
        ixs=[int(x) for x in xs]
        print("partitioning: Nquant,Nunit=",Nquant,Nunit)
        Ps=np.partition(allP,ixs)[ixs]
        print("done")
        quantiles=np.array([xs/NP,xs/Nunit,Ps])
        return np.array(quantiles)


# TODO, make this function also dependent on r_bar
def calcLogP(segments, populations, theta_p = [1, 1, 1], normalized = False, set_r_bar_fac=None):

        # Since we changed the number of populations we are considering,
        # This changed some data structure, this accounts for that change
        
        index = getIndex(populations)
        """
        if len(populations) == 3:
                index = [0, 1, 3]
        else:
                index = [0, 1, 2, 3]
        """

        #Allow r_bar to vary in the case that there is an extra value in theta_p
        if set_r_bar_fac is None:
                r_bar_factor=1.0
        else:
                r_bar_factor=set_r_bar_fac
        if len(theta_p) - len(populations) == 1:
                r_bar_factor=theta_p[-1]
                theta_p=theta_p[:-1]

        # TODO Check these constants
        T_alpha = 1648
        T_yr = 31557600.0
        T_total = T_alpha * len(segments)

        A_LPF = 10.37738
        K = T_alpha * A_LPF / T_yr #JGB: jacobian for impact dof accounted for in sums

        # pmax / pmin 
        # (Since both model and data are log binned, doesn't matter
        const = 1 #2 * (np.log(1000 / 0.1))

        norms = np.asarray([pop.norm for pop in populations]) #JGB: These are in units: all-sky per year
        #JGB What we really want here is the norm over the domain where the prior has support (and is constant)
        ### I haven't been able to determine the range of prior support by examining the MCMC code, but the
        ### cut-off in the range of the output values, as processed here is near 0.132.  This is not far from
        ### the minimum bin value of 0.1 in the population data, so may not be terrible to equivocate here.
        Fnorm = np.sum(np.multiply(norms, theta_p))
        if not normalized:
                # The total fluxes times our guess for theta_p 
                r_bar = Fnorm * K
                theta_p_norm = np.divide( theta_p, Fnorm )
        else:
                # detected / num_segments = # detected * Tobs / T_tot
                r_bar = N_DETECTED * T_alpha / T_total
                # Sophie: It doesnt matter that we are normalizing theta p, 
                # this is actually to normalize fluxes, I just want it 
                # outside the loop
                theta_p_norm = np.divide(theta_p, norms)  #JGB -changed name for clarity
        r_bar*=r_bar_factor

        logp = 0
        i = -1
        start_loop = time.time()
        #print('theta=',theta_p)
        #print('seg0 ifluxes=',segments[0].flux)
        #print('seg0  fluxes=',segments[0].flux[index, :])
        #print(segments[0].flux.shape)
        #print('->', np.matmul(theta_p, segments[0].flux[index, :]))

        for segment in segments:
                i += 1

                # Initialize segment lengths
                N_tot = segment.N
                N_1 = segment.Nimpact
                N_0 = N_tot - N_1
                sixteenoverlog10=16/np.log(10)
                if N_1 > 0:
                        # calc sum r_hat(psi_s)
                        r_seg_psi_theta = K * np.matmul(theta_p_norm, segment.flux[index, :])
                        #JGB: fluxes come in units per (2 deg)**2 per year
                        #JGB: For the momentum dimension the original prior was 1/8 dlnP, the new fluxes have units of either per x2 or per x5 even/odd
                        ### which we approximate as dlnp/(ln(10)/2), thus we have a lnP prior volume factor of 16/ln(10)
                        jacobian = 180**2/np.pi/np.cos(segment.lat[:]*np.pi/180.)  #JGB: provides all-sky normalizations for Ptot-integrated fluxes
                        jacobian *= sixteenoverlog10 # = 16/np.log(10) #lnP prior factor
                        r_hat = r_seg_psi_theta *jacobian
                else:
                        r_hat = [0]

                #print(i,'jac.mean,rh.min,rh.max,rh.mean',np.mean(jacobian),np.nanmin(r_hat),np.nanmax(r_hat),np.mean(r_hat))

                # ------ n == 1 ------
                n_1_term = r_bar * (1 / float(N_tot)) * np.sum(r_hat) / const
                n_0_term = (1 - (N_1 / float(N_tot))) * (1 - r_bar)
                
                if n_0_term + n_1_term <= 0:
                        print("Error! Prob is less than 0")
                #print(i,'n0term, n1term, r_bar, sum(r_hat)',n_0_term,n_1_term,r_bar,np.sum(r_hat))

                logp += np.log(n_0_term + n_1_term) # + prior + constant

        stop_loop = time.time()
        print(stop_loop - start_loop, 'seconds')

        print('Total probability ', logp)

        return logp

def writeLine(theta_p, prob, open_file):
        theta_p_string = ''
        for theta in theta_p:
                theta_p_string += '%e\t'%(theta)
        file_string = '%e\t'%(prob) + theta_p_string + '\n'

        open_file.write(file_string)

def mcmc(segments, populations, theta_p_init = [1, 1, 1, 1], N = 1000):
        mcmc_output = open(BASE_DIR + '/data/mcmc_output.txt', 'w+')

        print('starting MCMC')

        theta_p = np.asarray(theta_p_init)
        prob = calcLogP(segments, populations, theta_p)

        for n in range(N):
                print(100 * n / N)

                # Draw guess for theta_p
                variance = 1
                theta_p_proposed = np.abs(np.random.normal(0, variance, 4) + theta_p)
                prob_proposed = calcLogP(segments, populations, theta_p_proposed)

                # Calculate Hastings ratio
                H = min(1, np.exp(prob_proposed - prob))
                # Make random guess
                alpha = np.random.uniform(0, 1)
                
                # Check hastings ratio, if higher probability, always takes
                if H > alpha:
                        theta_p = theta_p_proposed
                        prob = prob_proposed

                # Write to file
                writeLine(theta_p, prob, mcmc_output) 
        mcmc_output.close()

        print('Done with MCMC')
        return

def makeGrid(N = 10, pop_num = 3):
        #Supports pop_num in [2,3]
        #JFC, HTC, AST, OCC

        # Cut up JFC first 
        # Initialize Grid of mixing probs
        Cs = np.ones((N ** (pop_num-1), pop_num))
        # v = (x + y) / 1
        # u = x / (x + y)

        # Help make uniform paramter space
        space = np.linspace(0, 1, N)
        Xs = []
        Ys = []
        Zs = []

        vspace=space
        if(pop_num<3):vspace=[1]
        for u in space:
                for v in vspace:
                        x = u * v / 1
                        y = v - x
                        z = 1 - (x + y)

                        Xs.append(x)
                        Ys.append(y)
                        Zs.append(z)

        #u = np.random.uniform(0, 1, N)
        #v = np.random.uniform(0, 1, N)
        """
        x = np.multiply(u, v)
        y = v - x
        z = 1 - (x + y)
        """

        Cs[:, 0] = Xs
        Cs[:, 1] = Ys
        if(pop_num==3):Cs[:, 2] = Zs

        return Cs

def getSegments(run_only_impacts = True, dataPath = BASE_DIR + dataDir,usePtot=True):

        df_veto = getVetoList()
        segments = []

        if run_only_impacts:
                index = df_veto.index[df_veto['isImpact']]
                impact_segs = df_veto['segment']
                only_impact_segments = df_veto['segment'][df_veto['isImpact']]
                for seg in only_impact_segments:
                        segments.append(simpleImpact(seg,usePtot=usePtot))

        else:
                filenames = []
                dirname=str(dataPath) + '/FLUX_IMPACTS_ALL'
                if not usePtot:dirname+='_INTP'
                for root, dirs, files in os.walk(dirname):
                        filenames = files
                all_segments = getGRSSegments(filenames)

                for seg in all_segments:
                        # Checks if impact is out of range and if the time 
                        # Contains an impact glitch
                        if isValid(seg, df_veto):
                                segments.append(simpleImpact(seg,usePtot=usePtot))
                        else:
                                continue

                # Checking lengths of segments
                #print("All segments :", len(all_segments))
                #print("Good segments:", len(segments))
                #good_segments = np.loadtxt(str(dataPath)+ '/segment_list.txt', dtype= int)
                #print("Jake segments", len(good_segments))
                

        return segments

def main(do_mcmc = False, save_filename = 'logp.txt', N = 100, normalized = True, usePtot=True,
                pop_names = ['JFC', 'HTC', 'AST', 'OCC', 'Uniform'], overwrite = False, run_only_impacts = True, set_r_bar_fac=None):
        
        # Set up directory structure
        dataPath = pathlib.Path(BASE_DIR + '/data')
        impact_dir = pathlib.Path(str(dataPath) + '/ALL_IMPACTS')
        #save_dir = pathlib.Path(str(dataPath) + '/FLUX_IMPACTS_ALL')
        modelDir = pathlib.Path(str(dataPath) + '/models')

        #usePtot = True
        save_file = str(dataPath) + '/' + save_filename
        
        if overwrite:
                write_file = open(save_file, 'w+')
        else:
                write_file = open(save_file, 'a+')

        start_readin = time.time()
        segments = getSegments(run_only_impacts = run_only_impacts, dataPath = BASE_DIR + dataDir,usePtot=usePtot)
        end_readin = time.time()
        print("Time to read segments in :", (end_readin - start_readin))
        dfrac = np.asarray([len(s.lon) / s.N for s in segments])
        dfrac = np.sort(dfrac)

        #generate quantiles
        quantiles=calcPQuantiles(segments)
        np.save('momentum_quantiles.npy', quantiles)

        fig, ax = plt.subplots()

        int_frac = dfrac
        int_frac[0] = dfrac[0]
        for i in range(len(dfrac) - 1):
                i += 1
                int_frac[i] = int_frac[i - 1] + dfrac[i]
        
                
        np.save('dfrac', dfrac)

        from copy import copy
        fig, ax = plt.subplots()
        dfrac = np.load('dfrac.npy')
        drfrac = np.sort(dfrac)

        int_frac = copy(dfrac)
        int_frac[0] = dfrac[0]
        for i in range(len(dfrac) - 1):
                i += 1
                int_frac[i] = int_frac[i - 1] + dfrac[i]

        ax.plot(dfrac, int_frac)
        ax.set_yscale('log')

        #return

        #Handle request to vary r_bar
        vary_r_bar=False
        if 'r_bar' in pop_names:
                pop_names.remove('r_bar')
                vary_r_bar=True
        
        # Read in populations
        print("Initializing Populations")
        pop_time = time.time()
        populations = []
        for p in pop_names:
                populations.append(pop(modelDir, p, usePtot))
        print("Time to read in populations", time.time() - pop_time)
        if vary_r_bar: pop_names+=['r_bar']

        print(1648 * len(segments), 'seconds observed')
        write_file = open(save_file, 'a+')
        write_file.write('logp\t')
        for o in pop_names:
                write_file.write('%s\t'%(o))
        write_file.write('\n')
        write_file.close()


        # MCMC would be good if we had a super computer
        if do_mcmc:
                #calcLogP(segments, populations, theta_p = [1, 1, 1, 1])

                theta_p_init = np.ones(len(populations))
                mcmc(segments, populations, theta_p_init = [1, 1, 1] , N = 100)

        else:
                theta_p_list = makeGrid(N = N, pop_num = len(populations))
                if vary_r_bar:
                        plist=theta_p_list.tolist()
                        theta_p_list=np.array([p+[r] for p in plist for r in np.arange(1,N+1)*4.0/N])
                count = -1
                for theta_p in theta_p_list:
                        count += 1
                        print(theta_p)
                        print(count * 100 / len(theta_p_list) , '% Complete')
                        logp = calcLogP(segments, populations, theta_p, normalized = normalized, set_r_bar_fac=set_r_bar_fac)
                        write_file = open(save_file, 'a+')
                        writeLine(theta_p, logp, write_file)
                        write_file.close()
        print("results in '"+save_file)
        return

#main(do_mcmc = False, save_filename = 'grid_smart_logp.txt', normalized = False)
#main(do_mcmc = False, save_filename = 'grid_normed.txt', normalized = True) #8 / 1 / 18
#main(do_mcmc = False, N = 100, save_filename = 'grid_normed_0802.txt', normalized = True) #8 / 2 / 18
#main(do_mcmc = False, N = 10, 
#               save_filename = '/population_ratios/grid_JFC_HTC_Uniform_0806.txt', normalized = True,
#               pop_names = ['JFC', 'HTC', 'Uniform'],
#               ) #8 / 6 / 18
#main(do_mcmc = False, N = 10,
#       save_filename = '/population_ratios/grid_OCC_JFC_HTC_0808.txt', normalized = True,
#       pop_names = ['OCC', 'JFC', 'HTC'],
#       ) #8 / 8 / 18


#Dec6:
#redefine(dataDir = '/data', saveDir = '/FLUX_IMPACTS_ALL')
###main(do_mcmc = False, N = 100, save_filename = 'population_ratios/grid_normed_20181206.txt', pop_names = ['JFC', 'HTC', 'Uniform'], normalized = True)
###main(do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_20181206.txt', normalized = True, pop_names = ['JFC', 'HTC', 'Uniform']) 
#main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_20190219.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform']) 
#main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_rbar2_20190219.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform'],set_r_bar_fac=2.0) 
#main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_rbar_20190219.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'r_bar']) 
#main(run_only_impacts=False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_OCC_JFC_HTC_20190219.txt', overwrite=True, normalized = True, pop_names = ['OCC', 'JFC', 'HTC']) 
#main(do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_OCC_JFC_HTC_20181221.txt', overwrite=True, normalized = True, pop_names = ['OCC', 'JFC', 'HTC']) 
###main(run_only_impacts = False, N = 10, save_filename = 'population_ratios/only_impacts_20181206.txt', normalized = True, pop_names = ['JFC', 'OCC', 'HTC'], overwrite = True)
##
#main(run_only_impacts=False, do_mcmc = False, N = 20, usePtot=False, save_filename = '/population_ratios/grid_intp_OCC_JFC_HTC_20190319.txt', overwrite=True, normalized = True, pop_names = ['OCC', 'JFC', 'HTC'])
#main(run_only_impacts = False, do_mcmc = False, N = 20, usePtot=False, save_filename = '/population_ratios/grid_intp_JFC_HTC_Uniform_20190319.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform'])
#main(run_only_impacts = False, do_mcmc = False, N = 20, usePtot=False, save_filename = '/population_ratios/grid_intp_JFC_HTC_rbar_20190319.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'r_bar'])
#March2019 usePtot-False test
#redefine(dataDir = '/data', saveDir = '/FLUX_IMPACTS_ALL_INTP',usePtot=False)
##
#rerun to verify code vs 190219 (start with rerun of redefine dir having saved off previous data to FLUX_IMPACTS_ALL.Dec)
redefine(dataDir = '/data', saveDir = '/FLUX_IMPACTS_ALL')
main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_20190327.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform']) 
main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_rbar2_20190327.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform'],set_r_bar_fac=2.0) 
main(run_only_impacts = False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_rbar_20190327.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'r_bar']) 
main(run_only_impacts=False, do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_OCC_JFC_HTC_20190327.txt', overwrite=True, normalized = True, pop_names = ['OCC', 'JFC', 'HTC']) 
#
#and also with rerun of
"""

good_segments = np.loadtxt(BASE_DIR + '/data/searchTimes0726.txt', dtype= int)
for s in good_segments:
        try:
                print(simpleImpact(s).segment)
        except:
                continue

        # Don't check segments I have vetoed
        if isValid(s, df_veto):
                segments.append(simpleImpact(s))
        else: 
                continue
        """




