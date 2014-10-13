import logging
from scipy.interpolate import interp1d
import math

version = '0.2.1'

# dictionary containing options
options = {}

def setup_logging(filename, write_on_screen, display_debug, logger_name):
    """Setup logging. Write logs to file with name specified in filename.
    If on_screen is true additionally write on screen.
    """
    # http://docs.python.org/library/logging.html
    logging.basicConfig(level = logging.DEBUG,
                        filename = filename,
                        filemode = 'w',
                        format = '%(levelname)-8s %(message)s')
    console_logger = logging.StreamHandler()

    if display_debug == True:
        console_logger.setLevel(logging.DEBUG)
    elif write_on_screen == True:
        console_logger.setLevel(logging.INFO)
    else:
        console_logger.setLevel(logging.WARNING)

    format_for_console = logging.Formatter('%(levelname)-8s %(message)s')
    console_logger.setFormatter(format_for_console)
    logging.getLogger(logger_name).addHandler(console_logger)


def parse_command_line_arguments(parser):
    arguments = parser.parse_args()
    args_dictionary = vars(arguments)

    # show version if asked to
    if args_dictionary['version'] == True:
        print 'Script version:', version
        exit(0)

    return args_dictionary


def check_option(option):
    return option in options.keys() and options[option] == 'True'


def transposed(matrix):
    """Return transposed matrix (list of lists).

    This function can handle non-square matrices.
    In this case it fills shorter list with None.

    >>> transposed( [[1,2,3], [3,4]] )
    [[1, 3], [2, 4], [3, None]]
    """
    return map(lambda *row: list(row), *matrix)


def interpolate(XYpairs, step = 0.01):
        
    xold = [x for (x,y) in XYpairs[0]]
    yold = [y for (x,y) in XYpairs[0]]
        
    spline = interp1d(xold, yold)
    
    xmax = max(xold)
    xnew = []
    xn = min(xold)
    while( xn < xmax ):
        xnew.append(xn)
        xn += step

    ynew = spline(xnew)
    
    return zip(xnew,ynew)


def make_function(XYpairs):
    """Return a function approximating data.

    This function takes list of (x,y) pairs and returns function of x.
    Between known xes linear interpolation is performed.

    TODO fix not working doctests, currently disabled
    ??? make_function( [(0,0), (1,2)] )(0.5)
    1

    The x values HAVE TO be evenly spread, if not error is thrown.

    ??? make_function( [(0,0), (1,2), (4,8)] )
    ERROR: Data in Bragg Peak database is not evenly spread. Please make it so.

    When trying to reach range greater then data 0 is returned.
    When trying to reach less than data last known is returned and
    warning message is displayed.

    ??? make_function( [(0,0), (1,2)] )(2)
    0
    ??? make_function( [(1,2), (2,3)] )(0)
    2
    """

    if len(XYpairs) == 0:
        logging.error('Trying to make_function() out of list of 0 length. Exiting.')
        exit(1)
    # TODO decide if function should do something reasonable when one pair is passed
    if len(XYpairs) == 1:
        logging.error('Trying to make_function() out of list of 1 length. Exiting.')
        exit(1)


    X = [ x for (x,y) in XYpairs]
    Y = [ y for (x,y) in XYpairs]

    # list of differences between subsequent Xes
    diffs = [X[i+1] - X[i] for i in range(len(X) - 1)]

    x_step = sum(diffs) / len(diffs)

    # TODO decide what to do if xes are the same
    if x_step == 0:
        logging.error('Argument to make_function has same xes in every pair. Exiting.')
        exit(1)


    # TODO:
    # replace following by something simpler, like 
    # diffs_reduced = [diff - diffs[0] for diff in diffs]
    # if ( max(diffs_reduced) == min(diffs_reduced) == 0 )...
    
    #check if x step is constant
    if len(XYpairs) == 2:
        standard_deviation_of_diffs = 0
    else:
        standard_deviation_of_diffs = ( sum([(d - x_step)**2 for d in diffs]) / (len(diffs) - 1) ) ** 0.5

    #TODO allow data not to be evenly spread
    if standard_deviation_of_diffs > 0.01 * x_step:
        logging.error('Data in Bragg Peak database is not evenly spread. Please make it so. Exiting.')
        exit(1)

    def func(x):
        if math.isnan(x):
            logging.warning('NaN appeared! AAAAAAAAA!')
            return float('nan')

        """function approximating data"""
        if x < X[0]:
            # it only dangerous when x > 0
            if x > 0:
                pass
                #logging.warning('Trying to get BP value from outside of range given in database. Given x: %f, Range: %f - %f. Returning last known value, i.e. %f', x, X[0], X[-1], Y[0])
            #return Y[0]
            return float('nan')
        elif x > X[-1]:
            return 0

        minX_i = int((x - X[0]) / x_step)
        if minX_i == len(X) - 1:
            return Y[minX_i]

        return Y[minX_i] + (x - X[minX_i]) * (Y[minX_i + 1] - Y[minX_i]) / x_step

    return func


def x_at_given_y(input_data, y_value, search_from_end = False, interpolate = True):
    """ TODO
    >>> x_at_given_y( [(0,0), (1,1), (3,3), (4,4)], 3.5 )
    3.5
    >>> x_at_given_y( [(0,0), (1,1), (3,3), (4,4), (5,3), (6,2), (8,0)], 1 , search_from_end = True)
    7.0
    >>> x_at_given_y( [(0,0), (1,1), (3,3), (4,4), (5,3), (6,2), (8,0)], 1 , search_from_end = True, interpolate = False)
    8.0
    >>> x_at_given_y( [(0,0), (1,2)], 2 )
    1.0
    >>> x_at_given_y( [(0,0), (1,2)], 0 )
    0.0
    """
    # sort comparing [0]th element from each pair
    input_data_sorted_by_x = sorted(input_data, key = lambda pair: pair[0])
    
    if search_from_end:
        input_data_sorted_by_x.reverse()
                
    if input_data_sorted_by_x[0][1] > y_value:
        return input_data_sorted_by_x[0][0]

    x_value = float('nan')
    for i in range(len(input_data_sorted_by_x)-1):
        xcur,ycur = input_data_sorted_by_x[i]
        xnext,ynext = input_data_sorted_by_x[i+1]
        if ynext == ycur and y_value == ycur:
            x_value = xcur
            break
        if y_value >= min(ycur,ynext) and y_value <= max(ycur,ynext):
            if interpolate:            
                x_value = xcur + ((xnext-xcur)/float(ynext-ycur))*(y_value-ycur)
            else:
                x_value = xcur
            break
        
    return float(x_value)


def calculate_BP_x_for_height_prox(BP_list, value_to_max):
    """Return position where BP first reaches point where
    it value to max ratio is as given.
    """
    # calculate maximum value
    values = transposed(BP_list)[1]
    maximum = max(values)

    x = x_at_given_y(BP_list, value_to_max * maximum, search_from_end=False)
    
    return x
    

def calculate_BP_x_for_height_dist(BP_list, value_to_max):
    """Return position where BP drops below point where
    it value to max ratio is as given.
    """
    # calculate maximum value
    values = transposed(BP_list)[1]
    maximum = max(values)

    x = x_at_given_y(BP_list, value_to_max * maximum, search_from_end=True)

    return x



def calculate_BP_80_width(BP_list):
    """Take list of (x,y) pairs representing a BP
    and return it's width at 80% of it's maximum.
    
    TODO doctest missing
    """
    # calculate maximum value
    values = transposed(BP_list)[1]
    maximum = max(values)
    value = 0.8 * maximum

    # first we go from the left
    left_x = x_at_given_y(BP_list, value, search_from_end=False)
    
    # now from the right
    right_x = x_at_given_y(BP_list, value, search_from_end=True)
        
    return right_x - left_x


def calculate_BP_range(BP):
    """
        TODO function, comment and doctest missing
    """

    return calculate_BP_x_for_height_dist(BP, 0.9)


def calculate_BP_maximum_position(BP):
    """Take list of (x,y) pairs representing a BP
    and return position of it's maximum.

    >>> calculate_BP_maximum_position([ (0,0), (1,1), (2,0) ])
    1
    """
    #sort pairs by y values
    sorted_XY = sorted(BP, key = lambda pair: pair[1])
    return sorted_XY[-1][0]

        

def calculate_BP_distal_falloff():
    pass


def calculate_BP_max_to_plateau():
    pass


def normalize_SOBP(SOBP_list):
    
    SOBP_begin = x_at_given_y(input_data = SOBP_list,
                                     y_value = 0.99,
                                     search_from_end = False)
    SOBP_end = x_at_given_y(input_data = SOBP_list,
                                   y_value = 0.99,
                                   search_from_end = True)

    if SOBP_end < SOBP_begin:
        return None
    
    plateau = [pair[1] for pair in SOBP_list if pair[0] > SOBP_begin and pair[0] < SOBP_end]
                
    normalizationFactor = sum(plateau)/len(plateau)
    
    SOBP_normalized_X = [ pair[0] for pair in SOBP_list]
    SOBP_normalized_Y = [ pair[1]/normalizationFactor for pair in SOBP_list]
    
    SOBP_normalized = zip( SOBP_normalized_X, SOBP_normalized_Y)
    
    return SOBP_normalized