#!/usr/bin/env python 
#!-*- coding: utf-8 -*-

'''Simplex Two-Phase Implementation'''

import numpy as np

__author__ = "Irshad Ahmad Bhat"
__version__ = "1.0"
__email__ = "irshad.bhat@research.iiit.ac.in"


class two_phase_simplex():

    '''Two Phase Simplex method operates on linear programs in standard form'''

    def __init__(self):

        self.table = None
        self.RHS = None
        self.var = None
        self.min_max = None
        self.constraint_count = None

    def get_equations(self):

        '''get equations and set initial values for parameters'''

        try:
            self.min_max, self.constraint_count = raw_input('\nInput:\n\n').split()
        except ValueError:
            print 'First line should contain two values'
            print 'exiting now...'
            exit()

        if self.min_max not in ['MIN', 'MAX']:
            print "No 'MIN'|'MAX' string found or not at 1st position"
            print "exiting now..."
            exit()

        try:
            self.constraint_count = int(self.constraint_count)
        except ValueError, e:
            print "Constraints count should be an integer, got '{}'".format(constraint_count)
            print 'exiting now...'
            exit()

        equations = list()
        RHS = list()

        # Get objective function
        eq = raw_input().split()
        self.var = len(eq)
        equations.append([int(i) for i in eq])

        # Get constraints
        for i in range(self.constraint_count):
            eq = raw_input().split()
            if len(eq) != (self.var + 1):
                print "coeffcient count mismatch"
                print "exiting now..."
                exit()
            try:
                RHS.append(int(eq[-1]))
                equations.append([int(i) for i in eq[:-1]])
            except ValueError, e:
                print e
                print "exiting now..."
                exit()

        # set RHS of objective function to 0
        RHS.append(0)
        equations.append(equations.pop(0))

        self.table = equations
        self.RHS = RHS

    def pprint(self, table):

        '''Pretty print a table to check intermediate results'''

        print
        for i in table:
            for j in i:
                print str(round(j,2))+'\t',
            print
        print

    def final_results(self, opt_point):

	'''post-processing of results'''

        # change sign of optimal value if MAX problem
        opt_value = self.table[-1][-1]
        if self.min_max == 'MAX':
            opt_value *= -1.0

        # get values of optimal point
        final_opt_pnt = list()
        for i in [0,1]:
            if i in opt_point:
                final_opt_pnt.append(self.table[opt_point.index(i)][-1])
            else:
                final_opt_pnt.append(0.0)

        return final_opt_pnt, opt_value

    def standard_form(self):

        '''
        Convert the given equations into their standard forms by adding 
        slack variables. Also change sign of optimization function if 
        it is a minimization problem so as to solve for maximization.
        '''

        # add slack varaibles
        for eq in self.table:
            eq.extend([0]*self.constraint_count)
        for i, eq in enumerate(self.table[:-1]):
            if self.RHS[i] < 0:
                eq[i+self.var] = -1
            else:
                eq[i+self.var] = 1

        # convert minimization problem to maximization by changing sign 
        self.table = np.array(self.table, dtype=float)
        if self.min_max == 'MIN':
            self.table[-1] = self.table[-1]*-1 + 0.0

        self.RHS =  np.array(self.RHS)
        self.table = np.hstack((self.table, self.RHS[:,np.newaxis]))

    def add_artf(self):

        '''Add artificial variables if there is no initial 
        basic feasible solution (BFS)'''

        bfs = True
        for i in range(self.var, self.var + self.constraint_count):
            if np.sum(self.table[:,i]) != 1:
                bfs = False
                art_var = np.zeros((self.constraint_count + 1, 1))
                art_var[i-self.var] = 1.0
                self.table = np.hstack((self.table, art_var))
                self.table[:,(-2,-1)] = self.table[:,(-1,-2)]

        return bfs

    def new_opt_fun(self):

        '''
        Create a new optimization function for the first phase and embed it 
        into the table. Also keep track of position of artificial variables
        '''

        artf_vars = list()
        for i in range(self.var + self.constraint_count, self.table.shape[1]-1):
            row_id = np.where(self.table[:,i]<0)[0]
            if list(row_id):    # NOTE np.array([0]) concidered None
                artf_vars.append(self.var + row_id[0])
                self.table[-1][self.var + row_id[0]] = -1
                self.table[-1] += self.table[row_id[0]]

        return artf_vars

    def simplex(self,opt_point):

        '''
        Loop until optimal point is reached or we get an infeasible solution
        '''

        self.pprint(self.table)
        while not np.all(self.table[-1][:-1] <= 0):
            pivot_col = np.argmax(self.table[-1][:-1])
            # check for feasibility
            if np.all(self.table[:,pivot_col][:-1] <= 0) or \
		np.any(self.table[:,pivot_col][:-1] == 0) or \
                np.all(self.table[:,-1][:-1] < 0):
                print 'Infeasible solution'
                exit()
            theta = self.table[:,-1][:-1] / self.table[:,pivot_col][:-1]
	    # check for feasibility
            if np.all(theta < 0):
                print 'Infeasible solution'
                exit()
            # set negative ratios to some large no. (+ infinity)
            theta[self.table[:,-1][:-1] < 0] = float('inf')
            theta[self.table[:,pivot_col][:-1] < 0] = float('inf')
            # get pivot and convert pivot row to 1 using row opeartion
            pivot_row = np.argmin(theta)
            self.table[pivot_row] /= self.table[pivot_row][pivot_col]
            # set Identity matrix using row operations
            for i in range(self.constraint_count+1):
                if i == pivot_row: continue
                self.table[i] = self.table[i]-(self.table[i][pivot_col]/ \
                        self.table[pivot_row][pivot_col])* \
                        self.table[pivot_row]
            self.pprint(self.table)
            # keep record of basic variables
            opt_point[pivot_row] = pivot_col

        return opt_point

    def two_phase(self):

        '''
        Check for initail BFS. If there is initial BFS call Simplex else 
        save original OPT function, add artificial variables, get new OPT 
        function and call first phase of simplex. If initial BFS obtained, 
        call simplex second phase with the original OPT function.
        '''

        # save original OP function
        Z = self.table[-1]
        bfs = self.add_artf()

        # initial BFS
        init_opt_pt = range(self.var,self.var + self.constraint_count)

        # if initial BFS
        if bfs:
            print "\nSimplex Method...."
            opt_point = self.simplex(init_opt_pt)
            return self.final_results(opt_point)

        # if no initial BFS
        # if any RHS is negative, change sign 
        neg_RHS = self.table[:,-1] < 0
        self.table[neg_RHS] = self.table[neg_RHS] * -1.0 + 0.0

        # set new OP function with artificial variables
        self.table[-1] = np.zeros(self.table.shape[1])
        artf_vars = self.new_opt_fun()

        # call first-phase
        print "\n1st Phase...."
        opt_point = self.simplex(init_opt_pt)

        # remove artificial variable set original OP function
        self.table = np.delete(self.table, artf_vars, axis=1)
        self.table[-1] = Z

        # check for identity matrix 
        for row, col in enumerate(opt_point):
            self.table[-1] -= self.table[-1][col] * self.table[row]

        # call 2nd-phase
        print "\n2nd Phase...."
        opt_point = self.simplex(opt_point)

        return self.final_results(opt_point)


if __name__ == '__main__':

    # initialize object
    simx = two_phase_simplex()

    # get equations and standarize them
    simx.get_equations()
    simx.standard_form()

    # call two-phase simplex method
    opt_point, opt_value = simx.two_phase()

    # print results
    print '\nOutput:\n'
    print 'Optimal point is {}'.format(opt_point)
    print 'Optimal value is {}'.format(opt_value)
