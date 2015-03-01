# Two Phase Simplex:

Two Phase Simplex tableau method for the linear programming model

## Input Format:

The first line of input will contain a string 'MAX' or 'MIN' followed by the number of constraints N. The next line will contain the coefficients for the objective function. The next N lines will contain the coefficients ai1 , ai2 , ...,  for the i th constraint.

## Constraint:

Note that this method assumes that all the constraints are of '<=' form. So, if a constraint is of '>=' form, convert it to '<=' form by changing the sign of the constraint. This doesn't affect the optimal point or the optimal value. For example:

    'x1 + 2x2 >= 14' is equivalent to '-x1 - 2x2 <= -14'

## How to use ??

    python simplex_2Ph.py

## Test Cases:

Below is a list of examples and the corresponding input to the code:

    ===================
    1. Maximize:	
    	20x1 + 10x2
        s.t:	
	x1  - x2 <= 1
	3x1 + x2 <= 7

	Input: 
	MAX 2
	20 10
	1 -1 1
	3 1 7

    ===================
    2. Maximize:    
	3x1 + 4x2
        s.t:	
	2x1 + x2  <= 600
	x1  + x2  <= 225
	5x1 + 4x2 <= 1000
	x1  + 2x2 >= 150

	Input:
	MAX 4
	3 4
	2 1 600
	1 1 225
	5 4 1000
	-1 -2 -150

    ===================
    3. Minimize:	
	x - x2
	s.t:	
	x1  + x2 >= 2
    	-x1 + x2 >= 1
	x2 <= 3

	Input:
    	MIN 3
	1 -2
	-1 -1 -2
	1 -1 -1
	0 1 3
    
    ===================
    4. Minimize:	
	2x1 + 3x2
	s.t:	
    	4x1 + 2x2 >= 12
	x1  + 4x2 >= 6

	Input:
    	MIN 2
	2 3
	-4 -2 -12
	-1 -4 -6

## Contact :
    Irshad Ahmad Bhat
    MS-CSE IIITH, Hyderabad
    bhatirshad127@gmail.com
    irshad.bhat@research.iiit.ac.in
