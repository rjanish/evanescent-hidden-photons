""" Process output from sample_effective_current """

import numpy as np


def read_effective_current_output(filenames, prefix=""):
    """ read output files of sample_effective_current """

    # param file format assumptions
    param_types = {"r_N"  :int, "r_start"  :float, "r_end"  :float,
                   "phi_N":int, "phi_start":float, "phi_end":float,
                   "z_N"  :int, "z_start"  :float, "z_end":float,
                   "m_N"  :int, "m_start"  :float, "m_end":float,
                   "atol":float, "rtol":float,
                   "length":float, "radius":float, "mode":str}
    num_params = len(param_types)
    num_header_lines = num_params + 4
        # assume param header occupies lines 1 to num_params, then
        # data identifer header occupies lines Nparam + 1 to num_params + 3,
        # so data begins on line Nparam + 4
    column_names = ["m", "r", "phi", "z",
                    "re_jr", "error_re_jr", "im_jr", "error_im_jr",
                    "re_jphi", "error_re_jphi", "im_jphi", "error_im_jphi",
                    "re_jz", "error_re_jz", "im_jz", "error_im_jz"]

    # read files
    if isinstance(filenames, str):
        filenames = [filenames] # handle single string filename input
    runs = {}
    for filename in filenames: # assumes a list of filenames
        label = filename.split(sep='.')[0][len(prefix):]
            # remove prefix and file ext
        runs[label] = {}
        # read params
        with open(filename) as file:
            for i in range(num_params):
                key, value = next(file).split()
                runs[label][key] = param_types[key](value)
                    # assume each line of param header contains a single
                    # param name and value separated by whitespae
        # read data
        runs[label]["rawdata"] = np.loadtxt(filename, skiprows=num_header_lines)
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
