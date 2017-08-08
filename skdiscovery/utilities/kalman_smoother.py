# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
# This software is part of the NSF DIBBS Project "An Infrastructure for
# Computer Aided Discovery in Geoscience" (PI: V. Pankratius) and 
# NASA AIST Project "Computer-Aided Discovery of Earth Surface 
# Deformation Phenomena" (PI: V. Pankratius)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# """@package kalman_smoother
# KalmanSmoother and utilities used to smooth gps time series data
# """

import numpy as np
import pandas as pd
import scipy.optimize as opt


def KalmanFilter(in_data, t, sigma_sq, R, Pinit, x0=0, invert=False, clipping=5):
    """
    Runs the kalman filter on data
    
    @param in_data Input data
    @param t: Correlation time
    @param sigma_sq: FOGM variance
    @param R: Noise variance
    @param Pinit: Initial variance
    @param x0: Intial updated state (default: 0)
    @param invert: Run the filter backwards (boolean flag)
    @param clipping: Clipping factor to use when computing cost functions

    @return the predicted state
    @return the predicted covariance
    @return the updated state
    @return the updated covariance
    @return C_hat, the sample innovation variance
    @return L, a different log variance cost function
    """

    data = np.copy(in_data)

    if invert == True:
        data = np.flipud(data)

    phi = np.exp((-1.0) / t)
    Q = sigma_sq*(1-np.exp(-2 / t))

    sample_size = len(data)

    updated_state = np.zeros(sample_size)
    updated_covariance = np.zeros(sample_size)
    
    predicted_state = np.zeros(sample_size)
    predicted_covariance = np.zeros(sample_size)

    updated_state[0] = x0
    predicted_state[0] = x0
    updated_covariance[0] = Pinit

    # If R is a scalar turn it into a series
    if np.isscalar(R):
        R = R * np.ones(sample_size)
    elif invert == True:
        R = np.copy(R)
        R = np.flipud(R)

    for i in range(1, sample_size):

        #prediction
        if invert == True:
            predicted_state[i] = updated_state[i-1] / phi
            predicted_covariance[i] =  (updated_covariance[i-1] + Q) / phi**2

        else:
            predicted_state[i] = phi * updated_state[i-1]
            predicted_covariance[i] = phi * updated_covariance[i-1]*phi + Q

        # update 
        if pd.isnull(data[i]):
            updated_state[i] = predicted_state[i]
            updated_covariance[i] = predicted_covariance[i]
            
        else:
            Kgain = predicted_covariance[i] / (predicted_covariance[i] + R[i])
            updated_state[i] = predicted_state[i] + Kgain * (data[i] - predicted_state[i])
            updated_covariance[i] = (1 - Kgain)*predicted_covariance[i]

    vk = data - predicted_state
    vk_std = np.nanstd(vk)
    clip_index = np.piecewise(vk, [~pd.isnull(vk),pd.isnull(vk)],
                              [lambda x: x > (clipping * vk_std), False]).astype(np.bool)

    vk[ clip_index ] = np.nan
    num_vk_nonnan = np.count_nonzero(~pd.isnull(vk))

    if t < 0 or sigma_sq < 0:
        Chat = np.nan
    else:
        Chat = (1 / num_vk_nonnan) * np.nansum(vk**2)

    Ck = predicted_covariance + R
    L = (1 / num_vk_nonnan) * (np.nansum(np.log(Ck)) + np.nansum(vk**2 / Ck))

    if invert == True:
        predicted_state = np.flipud(predicted_state)
        predicted_covariance = np.flipud(predicted_covariance)
        updated_state = np.flipud(updated_state)
        updated_covariance = np.flipud(updated_covariance)

    return(predicted_state, predicted_covariance, updated_state, updated_covariance, Chat, L)


def FitFOGMParameters(data, Pinit=100, R=1, method='brute', x0=0, clipping=5):
    """
    Find best FOGM parameters for a given data set

    @param data: input data
    @param Pinit: Initial updated covariance
    @param R: Noise Variance
    @param method: Method used to fit FOGM parameters. Use "simple", "brute", or "igrid". 
    @param x0: Initial value of x0 to use in the kalman filter
    @param clipping: Clipping factor used when computing cost functions

    @return best fit correlation time
    @return FOGM variance
    @return Noise variance
    @return correlation time from L
    @return FOGM variance from Chat
    """

    # Wrapper to return C_hat cost function
    def CostWrapperC(t_sigma_sq, data, R, P, x0=0, clipping=5):
        return(KalmanFilter(data, t_sigma_sq[0], t_sigma_sq[1], R, P, x0=x0, clipping=clipping)[4])

    # Wrapper to return L cost function
    def CostWrapperL(t_sigma_sq, data, R, P, x0=0, clipping=5):
         return(KalmanFilter(data, t_sigma_sq[0], t_sigma_sq[1], R, P, x0=x0, clipping=clipping)[5])

    # bounds for minimizing data
    bounds = [(1,len(data)), (1, 2*np.std(data)**2)]


    # Minimizing using default algorithm
    if method == 'simple':
        # initial estimate for minimum point 
        initial_estimate = (100, 50)

        (min_t_Chat, min_sigma_sq_Chat) = opt.minimize(CostWrapperC, initial_estimate, args=(data, R, Pinit,x0, clipping), bounds=bounds).x
        (min_t_L, min_sigma_sq_L) = opt.minimize(CostWrapperL, initial_estimate, args=(data, R, Pinit, x0, clipping), bounds=bounds).x
        
    # Minimizing by grid search
    elif method == 'brute':
        (min_t_Chat, min_sigma_sq_Chat) = opt.brute(CostWrapperC, bounds, Ns=40, args=(data, R, Pinit, x0, clipping))[0:2]
        (min_t_L, min_sigma_sq_L) = opt.brute(CostWrapperL, bounds, Ns=40, args=(data, R, Pinit, x0, clipping))[0:2]


    # Use iterative grid method described in
    # http://hdl.handle.net/1721.1/69466 (Kang Hyeun Ji 2011 MIT PhD Thesis)
    elif method == 'igrid':

        # Initial search space
        search_space = bounds

        # bounds used for resetting search session
        reset_bounds = [(search_space[0][0], None), (0, None)]

        def find_fogm_params(f, search_space, data, R, Pinit, x0, clipping, reset_bounds, tol):

            # center point of search space
            t_center = np.average(search_space[0])
            s_center = np.average(search_space[1])

            # maximum number of search resets
            max_search_iters = 10

            # Run Iterative Grid Search (IGS) once to get initial results
            (results, cost_minimum)  =  IterativeGridSearch(f, (data, R, Pinit, x0, clipping), search_space, bounds = reset_bounds, tol=0.1)
            t = results[0]
            sigma_sq = results[1]


            for search_iter in range(1, max_search_iters):

                # Stop search if either parameter is outside of parameter space
                # This check is not needed if updated search space is restricted
                # from extending outside of original search space
                if t > search_space[0][1]:
                    break;
                if sigma_sq > search_space[1][1]:
                    break

                # Previous IGS results
                prev_t = t
                prev_sigma_sq = sigma_sq
                prev_cost_minimum = cost_minimum


                # Determine half length of new search space so that it doesn't extend passed boundries
                t_side = (search_space[0][1] - t_center) -  abs(t - t_center)
                s_side = (search_space[1][1] - s_center) -  abs(sigma_sq - s_center)

                # Update search space
                new_search_space = [(t-t_side, t+t_side), (sigma_sq - s_side, sigma_sq + s_side) ]
                 

                # Rerun IGS
                (results, cost_minimum) =  IterativeGridSearch(f, (data, R, Pinit, x0, clipping),
                                                               new_search_space, bounds = reset_bounds,
                                                               prev_minimum = prev_cost_minimum, tol=tol)
                t = results[0]
                sigma_sq = results[1]


                # Stop search if the estimated FOGM parameters haven't changed between searches
                if(abs(t - prev_t) < 1 and abs(sigma_sq - prev_sigma_sq) < 0.1):
                    break

            # Make sure correlation time is within original bounds
            if t > search_space[0][1]:
                t = search_space[0][1]

            # Make sure sigma squared is within original bounds
            if sigma_sq > search_space[1][1]:
                sigma_sq = search_space[1][1]

            return(t, sigma_sq)

        # Run searches on both cost functions
        (min_t_Chat, min_sigma_sq_Chat) = find_fogm_params(CostWrapperC, search_space, data, R, Pinit, x0, clipping, reset_bounds, tol=0.1)
        (min_t_L, min_sigma_sq_L) = find_fogm_params(CostWrapperL, search_space, data, R, Pinit, x0, clipping, reset_bounds, tol=0.1)
                
    return(min_t_Chat, min_sigma_sq_L, R*(min_sigma_sq_L / min_sigma_sq_Chat), min_t_L, min_sigma_sq_Chat)


def IterativeGridSearch(f, args, intervals, max_iter=50, tol=0.1, bounds = None, prev_minimum = None, verbose = False):
    """
    Find the minimum of f using an iterative grid search with 3 points per dimension

    @param f: Function to be minimized. The function must accept a tuple with coordinates for the first input.
    @param args: additional arguments to pass on to the function.
    @param intervals: Space that contains the minimum. Must be a list of tuples, even if only 1 dimension.
    @param max_iter: Maximum number of iterations before stopping search.
    @param tol: Error tolerance on result.
    @param bounds: Additonal set of bounds for ending search.
    @param prev_minimum: Previous minimum of function. If the current minimum is close to the previous minimum the serach will stop
    @param verbose: Output debugging information.
    
    @return A tuple containing a numpy array with the location of the minimum; and the minimum value of the function.
    """

    # Given a function, arguments for the function, and set of intervals,
    # this will return the subinterval that contains the minimum
    def find_subinterval(f, args, intervals):
        # for each interval, determine location of grid points
        interval_bounds = []
        for interval in intervals:
            low = (1.5*interval[0] + 0.5*interval[1]) / 2
            high = (0.5*interval[0] + 1.5*interval[1]) / 2
            interval_bounds.append( (low, high) )

        # evaluate grid points using scipy's brute function
        grid_eval = None
        
        if args == None:
            grid_eval =  opt.brute(f, interval_bounds, Ns=3, finish=None, full_output=True)[3]
            
        else:
            grid_eval =  opt.brute(f, interval_bounds, args=args, Ns=3, finish=None, full_output=True)[3]

        if verbose:
            print('start: ')
            x_locations = np.asarray( [ interval_bounds[0][0], np.average(interval_bounds[0]), interval_bounds[0][1] ] )
            y_locations = np.asarray( [ interval_bounds[1][0], np.average(interval_bounds[1]), interval_bounds[1][1] ] )

            for i in range(0,3):
                for j in range(0,3):
                    print('(', x_locations[i], ', ', y_locations[j], ') : ', grid_eval[i,j], sep='')


            print(end='\n\n')
              
        # determine index of minimum value
        min_index = np.unravel_index(np.nanargmin(grid_eval), tuple((np.ones(len(intervals)) * 3).astype(int)))

        minimum_value = grid_eval[min_index]

        # function that returns appropriate subinterval given the index and interval
        def calc_return_interval(index, interval):
            low = (1.5*interval[0] + 0.5*interval[1]) / 2
            high = (0.5*interval[0] + 1.5*interval[1]) / 2
            center = (interval[0] + interval[1]) / 2
            
            if index == 0:
                return((interval[0], center))
            elif index == 1:
                return((low, high))
            elif index == 2:
                return((center, interval[1]))
 
        # return subintervals that contain the minimum
        return(np.asarray(
            [ calc_return_interval(min_index[i], intervals[i]) for i in range(0, len(intervals)) ]),
               minimum_value )

    # function to check when all intervals are below an absolute tolerance
    def below_tol(intervals, tol=0.1):
        for interval in intervals:
            if(abs(interval[1] - interval[0]) > tol):
                return False

        return True

    # Function that performs and additional bounds check
    def in_bounds(bounds, intervals):
        for i in range(0, len(intervals)):
            res = (intervals[i][0] + intervals[i][1])/2
            if bounds[i][0] != None and res < bounds[i][0]:
                return False
            if bounds[i][1] != None and res > bounds[i][1]:
                return False
        return True
            
    func_minimum = None
    subintervals = intervals

    # keep dividing the original intervals until max iterations or
    # tolorance threshold is reached
    for total_iter in range(0, max_iter):
        prev_subintervals = subintervals

        (subintervals, func_minimum) = find_subinterval(f, args, subintervals)

        # Check if result is within tolerance
        if below_tol(subintervals, tol):
            if verbose:
                print("below tolerance:", total_iter)
            break

        # Check if result is outside of additional bounds check
        if total_iter > 0 and bounds != None and in_bounds(bounds, subintervals) == False:
            subintervals = prev_subintervals
            if verbose:
                print("Result is out of bounds")
            break

        # Check if result is similar to previous tolerance
        if prev_minimum != None and abs(func_minimum - prev_minimum) < tol:
            if verbose:
                print("Reached result very close to previous minimum")
            break

    # return the center point for each sub interval
    return( np.asarray([ ((interval[0] + interval[1]) / 2) for interval in subintervals ]), func_minimum)


def KalmanSmoother(in_data, Pinit=1e6, Restimate=1, clipping=5, method='simple', t=None, sigma_sq=None, R=1, verbose=False, max_clip_iter=10):
    """
    Smoother based on a forward and a backward kalman filter

    @param in_data: Data to be smoothed (must be in a Pandas DataFrame)
    @param Pinit: Initial updated covariance
    @param Restimate: Initial estimate for noise variance
    @param clipping: Iteratively remove points beyond clipping * MSE.
    @param method: Method used to fit FOGM parameters. Use either "simple", "brute", or "igrid". 
    @param t: Fixed correlation time to use. Both sigma_sq and R must also be specified.
    @param sigma_sq: Fixed sigma squared to use. Both t and R must also be specified.
    @param R: Fixed measurement error to use Both t and sigma_sq must also be specified.
    @param verbose: Output additional information.
    @param max_clip_iter: Maximum number of clip iterations.

    @return values smoothed by the kalman smoother
    @return associated variance of smoothed result
    @return t, same as input, might have been altered by fitting parameters
    @return sigma_sq, same as input, might have been altered by fitting parameters
    @return R, same as input, might have been altered by fitting parameters
    """

    data = pd.DataFrame.copy(in_data)

    if t == None or sigma_sq == None:
        fit = True
    else:
        fit = False

    for clip_iter in range(0, max_clip_iter+1):

        # start_value = data[np.argmax(~pd.isnull(data))]
        start_value = np.nanmedian(data)

        # determining optimal parameters for kalman smoother
        if fit == True:
            if verbose:
                print("Fitting FOGM parameters")
            (t, sigma_sq, R) = FitFOGMParameters(data, Pinit, Restimate, R=R, method=method, x0=start_value, clipping=clipping)[0:3]
        else:
            if verbose:
                print('Using fixed parameters for FOGM fitting')
              
        if verbose:
            print(t, sigma_sq, R)

        # Run the kalman filter forward
        (forward_states, forward_covariance) = KalmanFilter(data, t, sigma_sq, R, Pinit, x0=start_value)[2:4]

        # Now run the filter backwards
        (backward_states, backward_covariance) = KalmanFilter(data, t, sigma_sq, R,
                                                              Pinit = forward_covariance[-1]*1E6,
                                                              x0 = forward_states[-1],
                                                              invert=True)[0:2]

        # Computer smoothed values from filters
        Kgain = forward_covariance / (backward_covariance + forward_covariance)
        xsmooth = pd.Series(forward_states + Kgain * (backward_states - forward_states), index = data.index)
        Psmooth = pd.Series(forward_covariance - Kgain*forward_covariance, index=data.index)

        # estimate total deviations squared between data and filter
        deviation = np.sqrt(np.nansum( (xsmooth-data)**2) / np.count_nonzero(~pd.isnull(data)))

        # Find all outliers
        outlier_index = data.index[ np.ndarray.flatten(np.asarray(np.abs(xsmooth - data)  > (deviation * clipping))) ]

        if clip_iter != max_clip_iter:
            if verbose:
                print("deviation:", deviation, len(outlier_index))
            if len(outlier_index) == 0:
                break
            else:
                data.loc[outlier_index] = np.nan

    return xsmooth, Psmooth, t, sigma_sq, R


def FOGM(size, t, sigma_sq, R):
    """
    Generates data from a First Order Gaussian-Markov process

    @param size: Number of data points
    @param t: Correlation time
    @param sigma_sq: FOGM variance
    @param R: Measurement variance

    @return Data generated from a FOGM
    """

    Q = sigma_sq * (1 - np.exp(-2 / t))

    phi = np.exp((-1.0)/t)    
    
    fogm_states = np.zeros(size)
    fogm_states[0] = np.random.normal(0, np.sqrt(Q))


    for i in range(1, size):
        fogm_states[i] = fogm_states[i-1]*phi + np.random.normal(0,np.sqrt(Q))

        
    measurements = fogm_states + np.random.normal(0, np.sqrt(R), size)
    
    return(measurements, fogm_states)

