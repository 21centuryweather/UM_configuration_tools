import pandas as pd
import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from pathlib import Path

def add_options(parser):
    """
    Process the command-line arguments
    """
    parser.description = 'Plot the memory usage of the Unified Model I/O server'

    parser.add_argument('-f','--input-file', 
        dest='input_file', 
        type=str, 
        help="The I/O input file, e.g. ioserver_log.00023")

    parser.add_argument('-m','--memory-limit',
        dest='memory', 
        type=int, 
        help="The I/O server memory limit in megabytes, i.e. the value of ios_buffer_size in [namelist:ioscntl]")
    

if __name__ == '__main__':

    parser = ArgumentParser()
    add_options(parser)
    
    if len(sys.argv) == 1:
        print (" INFO : No input data supplied.")
        print (" INFO : usage: python %prog [options].")
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    input_file = '~/cylc-run/u-dq126/work/20220226T0000Z/Flagship_ERA5to1km_1km_RAL3P2_um_fcst_000/ioserver_log.00023'

    INPUT_DIR=Path('/home/548/pag548/cylc-run/u-dq126/work/20220226T0000Z/Flagship_ERA5to1km_1km_RAL3P2_um_fcst_000/')

    log_files = INPUT_DIR.glob('ioserver_log.*')

    all_data = []

    for log_file in list(log_files):
        print (f' INFO : Opening {log_file.name}')

        data = pd.read_csv(log_file,
                        sep='\s+',
                        header=None,
                        names=['time',log_file.name,'change'])

        data.set_index('time',inplace=True)
        all_data.append(data[log_file.name])

    fig,ax = plt.subplots(1,1)
    for data in all_data:
        data.plot(ax=ax)

    # Filter only the largest ones
    for data in all_data:
        if data.max() > 1e10:
            data.plot(ax=ax)
            print (data.name)