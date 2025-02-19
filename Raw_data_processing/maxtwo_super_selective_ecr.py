#!/usr/bin/env python
# coding: utf-8
import h5py
import numpy as np
import pandas as pd
import pickle as pkl
import os
import sys
import glob
import time
import traceback
try:
    import ruamel.yaml as yaml
except:
    import ruamel_yaml as yaml
from datetime import datetime
# import matplotlib.pyplot as plt
# import matplotlib as mpl
from pprint import pprint
from oitools import ecr

filename_config = sys.argv[1]

print(f'Reading from configuration file: {filename_config}')
with open(filename_config, 'r') as f:
    yaml_loader = yaml.YAML(typ='safe', pure=True)
    config = yaml_loader.load(f)

# Path of source .h5 files
path_source_files = config['paths']['source_files']
if not path_source_files.endswith('/'):
    path_source_files += '/'

# Path where results will be stored
path_results = config['paths']['results']
if not path_results.endswith('/'):
    path_results += '/'

# Get a list of files that we'll analyze (h5 files in the source directory or its subdirectories)
filenames = glob.glob(f'{path_source_files}/**/*.h5', recursive=True)
n_chars = len(path_source_files)
filenames = [f[n_chars:] for f in filenames]
filenames.sort()
filenames.reverse()

print(f'Found {len(filenames)} h5 files to process.')

# Create an object that will run the ECR super-selective algorithm using the specified config parameters
ss_ecr = ecr.SuperSelective(epsilon=config['super_sel']['epsilon'],
                            T_list=config['super_sel']['T_list'],
                            sigma_list=config['super_sel']['sigma_list'],
                            amp_thresh_percentile=config['data']['corr_amp_thresh_percentile'],
                            amp_thresh_std=config['data']['corr_amp_thresh_std'],
                            n_corr_peaks_max=config['super_sel']['n_corr_peaks_max'],
                            raster_dur=config['super_sel']['raster_dur'],
                            corr_type=config['super_sel']['corr_type'],
                            adj_threshold=config['super_sel']['adj_threshold'],
                            time_resolution=1/config['data']['fs'],
                            use_multiprocessing=True,
                            verbose=True)

os.makedirs(path_results, exist_ok=True)
datetime_start = datetime.now().strftime("_%Y%m%d_%Hh%Mm")


# Analyze data for all files and all wells (all organoids)
for i_file, filename in enumerate(filenames):
    print(f'Processing file {i_file+1} of {len(filenames)}, {filename}...')

    filename_results = os.path.join(path_results, filename[:-3]+datetime_start+'.pkl')
    filename_results_dir, _ = os.path.split(filename_results)
    os.makedirs(filename_results_dir, exist_ok=True)
    
    if os.path.exists(filename_results) and not config['super_sel']['recompute']:
        # A results file already exists. Don't recompute.
        print(f"Results file already exists. Skipping {filename}.")
        continue
        
    data_to_save = {}
    data_to_save['source_filename'] = filename
    data_to_save['config'] = config
    
    fullname = os.path.join(path_source_files, filename)
    ephys_file = h5py.File(fullname, 'r')
    well_keys = list(ephys_file['recordings']['rec0000'].keys())

    for i_well, well in enumerate(well_keys):
        try:
            t_process_start = time.time()
            print(f'Processing well {i_well+1} of {len(well_keys)}, {well}...')
            df_spikes = pd.DataFrame(np.array(ephys_file['recordings']['rec0000'][well]['spikes']))

            # Confirm that the sampling rate matches what is given in the config file
            fs_h5 = ephys_file['recordings']['rec0000'][well]['settings']['sampling'][0]
            assert config['data']['fs']==fs_h5, 'Sampling rate of h5 file'

            # Determine spike amplitude threshold
            spike_amp_thresh = np.percentile(df_spikes['amplitude'], config['data']['spike_amp_thresh_percentile'])
            print(f"{config['data']['spike_amp_thresh_percentile']}% spike amplitude percentile is {spike_amp_thresh}.")

            # Filter spikes by amplitude
            df_spikes = df_spikes[df_spikes['amplitude'] > spike_amp_thresh]

            # Get spike times
            df_spikes['time_sec'] = df_spikes['frameno'] / config['data']['fs']
            t_sig_start = np.min(df_spikes['time_sec'])
            df_spikes['time_sec'] = df_spikes['time_sec'] - t_sig_start
            t_sig_last = np.max(df_spikes['time_sec'])
            print(f'Recording duration is {t_sig_last:0.4f} seconds.')

            # Aggregate spikes into separate channels
            chan_nums = np.unique(df_spikes['channel'])
            t_spikes_in_chan = []
            n_spikes_in_chan = np.zeros(len(chan_nums))
            for i_chan, chan in enumerate(chan_nums):
                t_spikes_in_chan.append(df_spikes[df_spikes['channel']==chan]['time_sec'].values)
                n_spikes_in_chan[i_chan] = len(t_spikes_in_chan[-1])

            # Put results in a dict and save as a pickle file
            data_to_save[well] = {}
            data_to_save[well]['spike_amp_thresh'] = spike_amp_thresh
            data_to_save[well]['channel_numbers'] = chan_nums
            data_to_save[well]['channel_spikes_per_sec'] = n_spikes_in_chan / t_sig_last

            ## Process ECR in overlapping windows
            if config['windows']['win_dur']=='None':
                t_win_start = [0]
                t_win_stop = [t_sig_last + 0.01]
                win_dur = t_sig_last
            else:
                win_dur = config['windows']['win_dur']
                overlap_dur = config['windows']['win_overlap_dur']
                t_win_start = np.arange(0, t_sig_last-overlap_dur, win_dur)
                t_win_stop = t_win_start + win_dur

            for i_win, (t_start, t_stop) in enumerate(zip(t_win_start, t_win_stop)):
                print(f'\nProcessing window {i_win+1} of {len(t_win_start)}.')
                t_spikes_in_chan_in_window = [[t for t in chan if t>=t_start and t<t_stop] for chan in t_spikes_in_chan]
                mean_spikes = np.mean([len(t) for t in t_spikes_in_chan_in_window])
                print(f'Mean spikes per second per channel is {mean_spikes/win_dur:0.2f}.')
                
                # Estimate Effective Connectivity Reconstruction
                adj_matrix_predicted, votes, corr_peaks = ss_ecr(t_spikes_in_chan_in_window)
    
                data_to_save[well][f'win_{i_win}'] = {}
                data_to_save[well][f'win_{i_win}']['adj_matrix_predicted'] = adj_matrix_predicted
                data_to_save[well][f'win_{i_win}']['votes'] = votes
                data_to_save[well][f'win_{i_win}']['corr_peaks'] = corr_peaks

            with open(filename_results, 'wb') as f:
                pkl.dump(data_to_save, f, protocol=pkl.HIGHEST_PROTOCOL)

            print(f"Done. Processing this well's data took {time.time()-t_process_start:0.2f} seconds.\n\n")
        except:
            print(f"Error while processing well {well} of fullname.")
            print(traceback.format_exc())
