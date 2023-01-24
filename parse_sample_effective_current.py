""" Process output from sample_effective_current """

import numpy as np


def parse_header_plus_array(filenames, param_types, prefix="", jump=0): 
    """
    Parse data file with header composed of key-value pairs, each on a single line
    and whitespace seperated, followed by a data array.  Input is a list of 
    filenames, a dict of header key names and associated datatype.  Output is 
    a dict of processed data, keyed first by filename (striped of prefix and extenstion)
    and then by param name and 'rawdata' containing the data array. jump is the
    number of rows after the end of the header at which the data array starts.  
    """
    num_params = len(param_types)
    num_header_lines = num_params + jump
        # assume param header occupies lines 1 to num_params, then some other
        # stuff occupies lines Nparam + 1 to num_params + jump - 1, and
        # the data begins on line Nparam + jump
    if isinstance(filenames, str):
        filenames = [filenames] # handle single string filename input
    runs = {}
    for filename in filenames: 
        label = filename.split('.')[0][len(prefix):] # remove prefix and file ext
        runs[label] = {}
        with open(filename) as file:
            for i in range(num_params):
                key, value = next(file).split()
                runs[label][key] = param_types[key](value)
        runs[label]["rawdata"] = np.loadtxt(filename, skiprows=num_header_lines)
    return runs


def read_effective_current_output(filenames, prefix=""):
    """ read output files of sample_effective_current """

    # param file format assumptions
    param_types = {"r_N"  :int, "r_start"  :float, "r_end"  :float,
                   "phi_N":int, "phi_start":float, "phi_end":float,
                   "z_N"  :int, "z_start"  :float, "z_end":float,
                   "m_N"  :int, "m_start"  :float, "m_end":float, 
                   "length":float, "radius":float, "omega":float, "mode":str,
                   "atol":float, "rtol":float,}
    column_names = ["m", "r", "phi", "z",
                    "re_jr", "error_re_jr", "im_jr", "error_im_jr",
                    "re_jphi", "error_re_jphi", "im_jphi", "error_im_jphi",
                    "re_jz", "error_re_jz", "im_jz", "error_im_jz"]

    runs = parse_header_plus_array(filenames, param_types, prefix, jump=4)
    for filename in filenames: # assumes a list of filenames
        label = filename.split('.')[0][len(prefix):] # remove prefix and file ext
        # process data
        shape_4d = (runs[label]["m_N"],
                    runs[label]["r_N"],
                    runs[label]["phi_N"],
                    runs[label]["z_N"])
                       # this is the order of data columns in the file
        for num, name in enumerate(column_names):
            runs[label][name] = (
                runs[label]["rawdata"][:, num].reshape(shape_4d))
                # put data in meshgrid format

        # construct other coordinates and values
        runs[label]["x"] = runs[label]["r"]*np.cos(runs[label]["phi"])
        runs[label]["y"] = runs[label]["r"]*np.sin(runs[label]["phi"])
        runs[label]["mag_jr"]   = np.sqrt(
            runs[label]["re_jr"]**2   + runs[label]["im_jr"]**2)
        runs[label]["mag_jphi"] = np.sqrt(
            runs[label]["re_jphi"]**2 + runs[label]["im_jphi"]**2)
        runs[label]["mag_jz"]   = np.sqrt(
            runs[label]["re_jz"]**2   + runs[label]["im_jz"]**2)
        runs[label]["jtotal"]   = np.sqrt(
            runs[label]["im_jr"]**2   + runs[label]["re_jr"]**2 +
            runs[label]["im_jphi"]**2 + runs[label]["re_jphi"]**2 +
            runs[label]["im_jz"]**2   + runs[label]["re_jz"]**2)
    return runs
    

def read_overlap_output(filenames, prefix=""):
    """ read output files of sample_overlap """
    param_types = {"Rs":float, "Ls":float, "mode_s":str, 
                   "Rd":float, "Ld":float, "mode_d":str, 
                   "omega":float, "sep":float, 
                   "m_start":float, "m_end":float,  "m_N":int, 
                   "atol":float, "rtol":float}
    return parse_header_plus_array(filenames, param_types, prefix, jump=4)