
"""
Gellispie Simulation of Cancer Growth
made for growth of the form: dn/dt = r*n^(1-alpha)+theta
where theta is a transition from a progentitor state    
@author: sophiemcgibbon
"""
from __future__ import division
import random
import math 
from numpy import genfromtxt
import numpy as np


# Simulation scale
sim_length = 1000
num_of_colonies = 100000
num_of_initial_cells = 1
max_N = 100000

# Rates of reaction
r0 = 1#replication rate
alpha = 0.7#exponent on crowding term
theta = 0.1#influx from progenitor cells

print "alpha, theta", alpha, theta

        
class gillespie:
    """
    Algorithm of a stochastic simulation.
    """
    def __init__(self, stem, cancer, r0, alpha, theta):
        """
        Description of varliables:
        alpha - crowding term
        r0 - replication rate
        theta - transition rate from progenitor state
        N_s - number of cells in the progenitor state
        N_c - number of replicating cancer cells
        N_tot - total number of cells in the simulation
        """
        self.alpha = alpha
        self.r0 = r0
        self.theta = theta
        self.N_s = sum(stem)
        self.N_c = sum(cancer)
        self.N_tot = sum(stem)+sum(cancer)

    def calc_totn(self, t, N_0, r0, alpha, theta):
        '''Calculate the mean feild result for the total number of cells
        Both power law growth dynamics and exponential growth dynamics are options
        '''
        #N_t = (alpha*t+(N_0**alpha))**(1./alpha)# sub exponential growth
        #N_t = N_0*np.exp(t)#exponential growth
        return(N_t)

    def probability(self, a0,a1, a, reaction, reaction_position):
        """
        Order of append commands for progenitors
        is important. Since code selecting the reaction
        relies on the order of append commands
        This calculates the arrays needed and does not accound for crowding
        Description of variables:
        a - list of reactions and their total liklihoods
        a1 - sum of all liklihoods for cancer cells
        a0 - sum of all liklihoods for progenitor cells (constant as they do not replicate or die)
        """
        colony_num = int(math.floor(reaction_position/2.0))
        if reaction == 0: #progenitor --> cancer
            a[reaction_position+1] += r0 #one more cancer cell
            a1+=r0
        elif reaction == 1: #cancer replacates
            a[reaction_position] += r0 #one more cancer cell
            a1 += r0
        else: #initialization of array, do nothing
            return a0, a1, a
        return a0, a1,a


    def choose_reaction(self, a0,a1, a, alpha, N_t):
        r = random.random()
        reaction_position = 0
        tot = 0.
        temp =a0+a1*(1./N_t**alpha)#total of probabilities for progentitors and for cancer cells, including crowding factor
        tempa = [x*(1./N_t**alpha) if i%2 else x for i,x in enumerate(a)]#modification of list a to include crowding factor
        # Grow total sum until position of a random arrow is found
        while tot < r*temp:
            tot += tempa[reaction_position]
            reaction_position += 1 #increase after to avoid calling out of bounds
            
        reaction_position -= 1  #decrease by one because had to increase after call in loop                   
        colony_num = int(math.floor(reaction_position/2.0))#find the colony number of the reaction
        reaction = reaction_position%2#find which reaction occured
        return colony_num, reaction, reaction_position

    def calculate_timestep(self, a0, a1, alpha, theta, N_t):
        """
        Calculate timestep sampled from exponential distribution.
        Since random.random() generates number in [0, 1)
        1.0 - random.random() will exclude 0 and generate (0,1]
        """
        r = 1.0 - random.random()
        temp =a0+a1*(1./N_t**alpha)
        return (1./temp)*math.log(1./r)

    def update_state(self, colony_num, reaction):
        '''
        Update the satate of the system given the previous reaction
        
        '''
        if reaction == 0:#progenitor --> cancer (alpha)
            cancer[colony_num] += 1
            self.N_c += 1.0
        elif reaction == 1:#cancer replicats (r0)
            cancer[colony_num] += 1
            self.N_c += 1.0
        
        self.N_tot = self.N_s + self.N_c
        return()
        
        
    def run(self, steps):
        """
        Perform stochastic simulation
        """
        #initailize reaction type as nothing
        reaction = -1
        reaction_position = 0
        
        # initialize total population
        self.N_tot = 0
        #initaiize probbilities
        a0 = 0
        a1 = 0
        a = []

        #Calculate probabilities
        for i in range(num_of_colonies):

            n = stem[i]
            a.append(n*theta/len(cancer))
            a0 += n*(theta/len(cancer))

            n = cancer[i]
            a.append(n*r0)
            a1 += n*r0

        print a0, a1
        name = str(alpha)
        #file to save results to
        myfile = open('theta_'+str(theta)+'alpha_'+str(alpha)+'b'+str(len(cancer))+'_r1.0.txt','w')
        t=0.000000000
        i=0
        n_nonzero1 = 0.
        barcodes = []
        self.N_s = sum(stem)
        self.N_c = sum(cancer)
        self.N_tot = self.N_s + self.N_c
        for j in range(len(stem)):
            barcodes.append(stem[j]+cancer[j])
        N_0 = sum(cancer)
        N_t = 1#initialize as 1 to avoid out of bounds call, is actually 0

        #Do one round of the simulation
        a0,a1, a = self.probability(a0,a1, a, reaction, reaction_position)
        colony_num, reaction, reaction_position = self.choose_reaction(a0,a1, a, alpha, N_t)
        r = 1.0 - random.random()
        timestep = (1./(a0+a1))*math.log(1./r)
        self.update_state(colony_num, reaction)
        t += timestep
        self.N_s = sum(stem)
        self.N_c = sum(cancer)
        self.N_tot = self.N_s + self.N_c
        count = 0
        timearr = [0]
        trajec = list(cancer)
        for i in range(sim_length):
                count +=1
                a0,a1, a = self.probability(a0,a1, a, reaction, reaction_position)
                N_t = self.N_c
                colony_num, reaction, reaction_position = self.choose_reaction(a0,a1, a, alpha, N_t)
                timestep = self.calculate_timestep(a0,a1,alpha, theta, N_t)
                self.update_state(colony_num, reaction)
                t += timestep
                
                #use "stack" to put extra column on list, append the time, then use time list as column and stacked list as data for each run
                if np.log2(count)%1 == 0:
                    timearr.append(t)
                    trajec = np.column_stack((trajec,cancer))
                    print 'log', np.log2(count), 'count', count, 't', t, 'N(t)', N_t, 'N_c', self.N_c, 'nonzero', np.count_nonzero(cancer)
                    time_str = str(int(t*100000)/100000.0)                   
                    myfile.write("%5.9f\t%8.2f\n" %(t, self.N_c))
                self.N_s = sum(stem)
                self.N_c = sum(cancer)
                self.N_tot = self.N_s + self.N_c
        df = pd.DataFrame(trajec, columns = timearr)
        df.to_csv('alpha'+str(alpha)+'_theta'+str(theta)+'b_'+str(len(cancer))+'.csv', index = False)
        df.to_csv('alpha'+str(alpha)+'_theta'+str(theta)+'_b'+str(len(cancer))+'.csv', index = False)
        df.to_csv('alpha'+str(alpha)+'_noH_theta'+str(theta)+'_b'+str(len(cancer))+'.csv', header = False, index = False)
        np.savetxt("alpha"+str(alpha)+'_theta'+str(theta)+'_b'+str(len(cancer))+".csv", cancer, delimiter=',')
        myfile.close()
                        
#Initiate colonies
stem = []
cancer = []
for num in range(num_of_colonies):
    stem.append(1)
    cancer.append(0)

sim = gillespie(stem, cancer, r0, alpha, theta)
sim.run(sim_length)
