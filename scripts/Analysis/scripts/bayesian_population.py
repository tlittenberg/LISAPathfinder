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


def redefine(dataDir = '/data', saveDir = '/FLUX_IMPACTS_ALL'):
        """
        Creates .npy files used in the analysis

        """
        dataPath = pathlib.Path(BASE_DIR + dataDir)
        usePtot = True
        
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

                        if impact.lon_sun is None:
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
                                        flux_JFC = populations[0].getFlux(lons, lats, Ptots, norm = False)
                                        flux_HTC = populations[1].getFlux(lons, lats, Ptots,  norm = False)
                                        flux_AST = populations[2].getFlux(lons, lats, Ptots, norm = False)
                                        flux_OCC = populations[3].getFlux(lons, lats, Ptots,  norm = False)

                                        # This number is meaningless with no norm but we dont have it for simplification reading in 
                                        flux_Uniform = populations[4].getFlux(lons, lats, Ptots,  norm = False) 
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
                        
                        print(impact.segment,lons[0:3],lats[0:3],flux_JFC[0:3],flux_HTC[0:3])
                        N = impact.N * np.ones(len(lons))
                        segment = impact.segment * np.ones(len(lons))

                        data_array = [impact.lon_sun, impact.lat_sun, impact.Ptot, 
                                                  flux_JFC, flux_HTC, flux_AST, flux_OCC, flux_Uniform, N, segment]

                        np.save(file_numpy, data_array)
                        write_file.write('saved' + '\n')
                        write_file.write('\n')

class simpleImpact:
        def __init__(self, segment, direc = '/data/FLUX_IMPACTS_ALL'):
                segment_direc = BASE_DIR + direc
                data = np.load(segment_direc + '/' + str(segment) + '_grs1.npy')
                print("reading",segment)
                print(data.shape)
                try:
                        self.lon = data[0, :]
                        self.lat = data[1, :]
                        self.Ptot = data[2, :]
                        # JFC, HTC, AST, OCC, Uniform
                        self.flux = data[3:8, :]

                        self.N = int(data[8, 0])
                        self.segment = int(data[9, 0])

                except:
                        self.lon = [None]
                        self.lat = [None]
                        self.Ptot = [None]
                        self.flux = [None]
                        self.N = int(data[8])
                        self.segment = int(data[9])

                print('flux[:,0:3]=',self.flux[:,0:3])
                print('flux[:,1000:1003]=',self.flux[:,1000:1003])
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


# TODO, make this function also dependent on r_bar
def calcLogP(segments, populations, theta_p = [1, 1, 1], normalized = False):

        # Since we changed the number of populations we are considering,
        # This changed some data structure, this accounts for that change
        
        index = getIndex(populations)
        """
        if len(populations) == 3:
                index = [0, 1, 3]
        else:
                index = [0, 1, 2, 3]
        """

        # TODO Check these constants
        T_alpha = 1648
        T_total = T_alpha * len(segments)

        A_LPF = 10.37738
        jacobian = 1 / (4 * np.pi) #(From sky distribution)
        K = T_alpha * A_LPF * jacobian

        # pmax / pmin 
        # (Since both model and data are log binned, doesn't matter
        const = 1 #2 * (np.log(1000 / 0.1))

        norms = np.asarray([pop.norm for pop in populations])
        if not normalized:
                # The total fluxes times our guess for theta_p 
                r_bar = np.sum(np.multiply(norms, theta_p))
        else:
                # detected / num_segments = # detected * Tobs / T_tot
                r_bar = N_DETECTED * T_alpha / T_total
                # It doesnt matter that we are normalizing theta p, 
                # this is actually to normalize fluxes, I just want it 
                # outside the loop
                theta_p = np.divide(theta_p, norms)

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
                if segment.lon[0] is None:
                        N_1 = 0
                else:
                        N_1 = len(segment.lon)
                N_0 = N_tot - N_1


                if N_1 > 1:
                        # calc sum r_hat(psi_s)
                        r_seg_psi_theta = K * np.matmul(theta_p, segment.flux[index, :])
                        
                        #r_seg_psi_theta = np.asarray([K * np.sum(np.multiply(theta_p, segment.flux[index, i])) for i in range(N_1)])
                        
                        # events seen / events expected
                        r_hat = r_seg_psi_theta / r_bar
                else:
                        r_hat = [0]
                
                # ------ n == 1 ------
                n_1_term = r_bar * (1 / float(N_tot)) * np.sum(r_hat) / const
                n_0_term = (1 - (N_1 / float(N_tot))) * (1 - r_bar)
                
                if n_0_term + n_1_term <= 0:
                        print("Error! Prob is less than 0")

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
        #JFC, HTC, AST, OCC

        # Cut up JFC first 
        # Initialize Grid of mixing probs
        Cs = np.ones((N ** 2, pop_num))
        # v = (x + y) / 1
        # u = x / (x + y)

        # Help make uniform paramter space
        space = np.linspace(0, 1, N)
        Xs = []
        Ys = []
        Zs = []
        for u in space:
                for v in space:
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
        Cs[:, 2] = Zs



        return Cs

def getSegments(run_only_impacts = True, dataPath = BASE_DIR + dataDir):

        df_veto = getVetoList()
        segments = []

        if run_only_impacts:
                index = df_veto.index[df_veto['isImpact']]
                impact_segs = df_veto['segment']
                only_impact_segments = df_veto['segment'][df_veto['isImpact']]
                for seg in only_impact_segments:
                        segments.append(simpleImpact(seg))

        else:
                filenames = []
                for root, dirs, files in os.walk(str(dataPath) + '/FLUX_IMPACTS_ALL'):
                        filenames = files
                all_segments = getGRSSegments(filenames)

                for seg in all_segments:
                        # Checks if impact is out of range and if the time 
                        # Contains an impact glitch
                        if isValid(seg, df_veto):
                                segments.append(simpleImpact(seg))
                        else:
                                continue

                # Checking lengths of segments
                #print("All segments :", len(all_segments))
                #print("Good segments:", len(segments))
                #good_segments = np.loadtxt(str(dataPath)+ '/segment_list.txt', dtype= int)
                #print("Jake segments", len(good_segments))
                

        return segments

def main(do_mcmc = False, save_filename = 'logp.txt', N = 100, normalized = True, 
                pop_names = ['JFC', 'HTC', 'AST', 'OCC', 'Uniform'], overwrite = False, run_only_impacts = True):
        
        # Set up directory structure
        dataPath = pathlib.Path(BASE_DIR + '/data')
        impact_dir = pathlib.Path(str(dataPath) + '/ALL_IMPACTS')
        save_dir = pathlib.Path(str(dataPath) + '/FLUX_IMPACTS_ALL')
        modelDir = pathlib.Path(str(dataPath) + '/models')

        usePtot = True
        save_file = str(dataPath) + '/' + save_filename
        
        if overwrite:
                write_file = open(save_file, 'w+')
        else:
                write_file = open(save_file, 'a+')

        start_readin = time.time()
        segments = getSegments(run_only_impacts = run_only_impacts, dataPath = BASE_DIR + dataDir)
        end_readin = time.time()
        print("Time to read segments in :", (end_readin - start_readin))
        dfrac = np.asarray([len(s.lon) / s.N for s in segments])
        dfrac = np.sort(dfrac)

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

        # Read in populations
        print("Initializing Populations")
        pop_time = time.time()
        populations = []
        for p in pop_names:
                populations.append(pop(modelDir, p, usePtot))
        print("Time to read in populations", time.time() - pop_time)


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
                count = -1
                for theta_p in theta_p_list:
                        count += 1
                        print(theta_p)
                        print(count * 100 / len(theta_p_list) , '% Complete')
                        logp = calcLogP(segments, populations, theta_p, normalized = normalized)
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
main(do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_JFC_HTC_Uniform_20181221.txt', overwrite=True, normalized = True, pop_names = ['JFC', 'HTC', 'Uniform']) 
#main(do_mcmc = False, N = 20, save_filename = '/population_ratios/grid_OCC_JFC_HTC_20181221.txt', overwrite=True, normalized = True, pop_names = ['OCC', 'JFC', 'HTC']) 
###main(run_only_impacts = False, N = 10, save_filename = 'population_ratios/only_impacts_20181206.txt', normalized = True, pop_names = ['JFC', 'OCC', 'HTC'], overwrite = True)

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




