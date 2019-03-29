from impactClass import impactClass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import pathlib
import pickle
import pandas as pd
import os
from mpl_toolkits.mplot3d import Axes3D
import scipy as scp


class population:
        """
                Class to store data related to population models in JFC and HTC_30um files

                TODO:
                - Make work with varying estimation for momentum
                - Add Probability function
        """
        ## Private Functions ## 
        def __init__(self, modelDir, pop_type, usePtot = True):
                """
                        modelDir = directory where file is located
                        pop_type = String, Population type, one of 'JFC', 'HTC', 'Uniform', 'AST', 'OCC'
                        usePtot = Whether to integrate momentum out
                """

                # Setting wheter to interpolate with Ptot
                self.usePtot = usePtot
                self.pop_type = pop_type

                if pop_type == 'Uniform':
                        dataFile = None
                elif pop_type in ['AST', 'JFC', 'HTC', 'OCC']:
                        dataFile = pop_type + '_impulse.dat'    
                else:
                        print('pop_type not valid:', pop_type)
                        print('Please use one of : \n\t-JFC\n \n\t-HTC\n \n\t-AST \n\t-OCC \n\t-Uniform\n')
                        raise(ValueError)

                self.grid = self.getGrid()
                
                # For uniform just take the Given Grid
                if pop_type == 'Uniform':
                        self.df = self.getGrid()

                        # Set flux_map
                        self.setFluxMap()
                        if usePtot:
                                self.norm = np.sum(self.getFlux(self.df['lon'].values, self.df['lat'].values, self.df['Ptot'].values, norm = False))
                        else:
                                self.norm = np.sum(self.getFlux(self.df['lon'].values, self.df['lat'].values, norm = False))
                        print('uniform norm=',self.norm)
                        return

                # Read in data and convert to dataframe
                self.df = self.__readSkymapData__(modelDir, dataFile)

                # Initialize Norm
                self.norm = 1

                # Set flux_map
                self.setFluxMap()

                # Set norm
                try:
                        self.norm = np.sum(self.getFlux(self.grid['lon'].values, self.grid['lat'].values, self.grid['Ptot'].values, norm = False))
                except KeyError:
                        print('Ptot has been deleted')
                        self.norm = np.sum(self.getFlux(self.grid['lon'].values, self.grid['lat'].values, norm = False))

                return

        def getGrid(self):
                """ Returns our Grid """
                grid_lon = np.arange(-179, 181, 2)
                grid_lat = np.arange(-89, 91, 2)
                grid_Ptot = np.asarray([1e-7, 5e-7, 1e-6, 5e-6, 
                                                                1e-5, 5e-5, 1e-4, 5e-4, 1e-3]) * 10 ** 6 # Micro Ns

                lon = []
                lat = []
                Ptot = []
                for i in grid_lon:
                        for j in grid_lat:
                                # Includes Momentum
                                if self.usePtot:
                                        for p in grid_Ptot:
                                                lon.append(i)
                                                lat.append(j)
                                                Ptot.append(p)
                                else:
                                        lon.append(i)
                                        lat.append(j)
                flux = np.zeros(len(lon))
                
                if self.usePtot:
                        grid = {'lon' : lon, 'lat' : lat, 'Ptot' : Ptot, 'flux' : flux}
                else:
                        grid = {'lon' : lon, 'lat' : lat, 'flux' : flux}

                df_grid = pd.DataFrame(grid)

                return df_grid

        def __readSkymapData__(self, modelDir, dataFile = "JFC_30um"):
                """
                        modelDir = directory where file is located
                        datafile = filename (JFC_30um or HTC_30um)
                """
                # modules
                import numpy as np

                # Creates dataframe with Ptot in it
                df = pd.read_csv(str(modelDir) + '/' + dataFile, 
                                                names = ['lon', 'lat', 'Ptot', 'flux'], sep = '\s+')
                
                df['rads'] = np.deg2rad(df['lat'])

                # Converts to MicroNs
                df.eval('Ptot = Ptot * (10 ** 6)', inplace=True)

                # DO not divide by cos, gets even flux for each term
                #df.eval('flux = flux * 1 / cos(rads)', inplace=True)

                del df['rads']

                return df

        def integratePtot(self, df = None):
                if df is None:
                        df = self.df

                try:
                        del df['Ptot']
                except KeyError:
                        print('Ptot has already been deleted!')
                # All have the same lon, lat so dont have to 1/cos(lat) for flux
                grouped = df.groupby(['lon', 'lat'])['flux'].sum().reset_index()
                df = grouped

                return df

        def getP_n1(self):
                ''' 
                Gets probability p( n = 1 | theta_p)
                ie) Integrate over grid
                '''

                return np.sum(self.df['flux'].values )
                                


        def interpolate(self, integrate_mom = False):

                # Integrated map
                df_int = self.integratePtot(self.df.copy())
                self.flux_map_int = scp.interpolate.NearestNDInterpolator((df_int['lon'].values,
                                                                                                                          df_int['lat'].values),
                                                                                                                          df_int['flux'].values)
                if self.usePtot and not integrate_mom:
                        print("Will interpolate on grid:")
                        print("(lon,lat):",(df_int['lon'].values,df_int['lat'].values))
                        try:

                                self.flux_map = scp.interpolate.NearestNDInterpolator((self.df['lon'].values,
                                                                                                                  self.df['lat'].values,
                                                                                                                  self.df['Ptot'].values),
                                                                                                                  self.df['flux'].values)

                        except KeyError:
                                print("Could not interpolate with Momentum")
                                print("Ptot has been integrated out")
                        else:
                                return


                if not self.usePtot:
                        self.flux_map = self.flux_map_int

                return

        def setFluxMap(self):
                """
                        Sets the flux map by interpolating data
                """
                if self.pop_type == 'Uniform':
                        # Create a power law for momentum
                        def power_func(x, a = 1e-07, b = -1):
                                # a and b are 
                                # average fit for flux vs momentum
                                return a * (x ** b)

                        df_grid = self.getGrid()
                        self.df = df_grid.copy()

                        # Creates a uniform flux
                        if self.usePtot:
                                self.df['flux'] = (power_func(self.df['Ptot'].values) / len(self.df.index))
                        else:
                                self.df['flux'] = 1.0;
                                
                        # Interpolate uniform flux map
                        self.interpolate()
                else:
                        self.interpolate()
                return

        def getFlux(self, lon, lat, Ptot = None, norm = True):
                """
                Returns nearest flux value for the set of longitudes and latitudes
                lon - scalar or arraylike, set of longitudes
                lat - scalar or arraylike, set of lattitudes
                Ptot -> only used if self.usePtot, set of Momenta

                """
                if Ptot is None:
                        # Uses map not interpolated with momentum
                        if norm:
                                return self.flux_map_int(lon, lat) / self.norm
                        else: 
                                return self.flux_map_int(lon, lat)  

                # Used map interpolated with momentum
                if self.usePtot:
                        if norm:
                                return self.flux_map(lon, lat, Ptot) / self.norm
                        else:
                                return self.flux_map(lon, lat, Ptot)

        def calc_like_impact(self, impact, norm = True):
                """
                Calculates likelihood that impact came from 1 population

                        pop = populationClass instance
                        impact = impactClass instance
                        
                        norm says whether or not to use normed probabiliy
                        normed probability will not give us a rate

                returns likelihood of impact (scalar)
                """


                # Get impact in right coordinates
                #impact = impact.findSkyAngles()
                #impact = impact.SCtoSun()
                #impact = impact.SuntoMicro()

                # Gets lon lat in correct frame
                lons = impact.lat_sun
                lats = impact.lon_sun
                Ptots = impact.Ptot

                if Ptots is None:
                        print('None')

                try:
                        N_1 = len(lons)
                except:
                        N_1 = 0

                # ------- For n == 1 -------
                # Calculates f_pop = p(n = 1, psi | theta_p) 
                N_tot = impact.N
                p_hat_1 = 0.5

                # Integral over parameters 
                if lons is None:
                        # Very small number, can't be zero
                        f_pop = 1e-218
                elif self.usePtot:
                        f_pop = self.getFlux(lons, lats, Ptots, norm = norm) / p_hat_1
                else:
                        f_pop = self.getFlux(lons, lats, norm = norm) / p_hat_1
                

                # ----- For n == 0 -------
                p_hat_0 = 0.5
                N_0 = N_tot - N_1

                # p(n = 0 | theta_p)
                # p(n = 0 | theta_p) = 1 - p(n=1|theta_p) 
                # Sum over all space
                f_pop_1 = np.sum(self.getFlux(self.grid['lon'].values, self.grid['lat'].values, 
                                                 Ptot = self.grid['Ptot'].values, norm = norm))
                f_pop_0 = (1 - f_pop_1) / p_hat_0

                # Sum p(n = 1, psi | theta_p)
                sum_pop = np.sum(f_pop)
                
                if N_1 == 0:
                        likelihood = (N_0 / N_tot) * f_pop_0
                else:   
                        likelihood = (N_0 / N_tot) * f_pop_0 + sum_pop / N_1

                return likelihood

                

