#!/usr/bin/env python3

import os
import raw_class
import time
import subprocess
import argparse

ZSTD_PATH = 'zstd'  # '~/bin/zstd'
data_dir = '/data/'

DEBUG = True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Olaf2Raw', description="This convert a waveolaf file in UDP packets simulating real-time.",
                                     formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=160))
    parser.add_argument('-dst_add', dest='dst_add', default="127.0.0.1",
                        help="destination addresse default is -dst_add=\"127.0.0.1\"")
    parser.add_argument('-dst_port', dest='dst_port', default=5586,
                        help="destination port default is -dst_port=5586")
    parser.add_argument('-duration', dest='duration', default=None,
                        help="UDP socket duration default is the file duration in real time")
    parser.add_argument('-udp_speed', dest='udp_speed', default=1,
                        help="UDP socket speed default 1.0 corresponding to real-time")
    parser.add_argument('INPUT_ARCHIVE', nargs='+', help="Path to the Archives")

    args = parser.parse_args()

    # wavfile (waveolaf format)
    # file = '/data/B2217+47_1242.undysputedbk1.2022-03-04T12:03:10.000.rawolaf_tmp'
    file = args.INPUT_ARCHIVE[0]
    # file = '/data/B2217+47_1242.undysputedbk1.2022-03-04T12:03:10.000.zst'
    # zstd -d /data/B2217+47_1242.undysputedbk1.2022-03-04T12:03:10.000.zst -o /data/B2217+47_1242.undysputedbk1.2022-03-04T12:03:10.000.rawolaf_tmp

    if (file.split('.')[-1] == 'zst'):
        uncompressed_file = data_dir + '.'.join(os.path.basename(file).split('.')[:-1]) + '.rawolaf_tmp'
        print("Decompression of %s to %s" % (file, uncompressed_file))
        cmd = ZSTD_PATH + ' -d ' + file + ' -o ' + uncompressed_file
        completed = subprocess.run(cmd, shell=True)
        file = uncompressed_file

    # initialisation of the raw object containing the methodes
    raw_object = raw_class.Raw([file],
                               verbose=True)
    file_start = raw_object.get_starttime()   # mjd
    file_stop = raw_object.get_stoptime()  # mjd
    file_duration = file_stop.unix - file_start.unix

    print(" File start at smjd : %d sec" % (file_start.unix % (24 * 3600)))
    print(" File stop  at smjd : %d sec" % (file_stop.unix % (24 * 3600)))
    print(" File duration is   : %d sec" % (file_duration))

    raw_object.start(speed=float(args.udp_speed),
                     adresse=(args.dst_add, int(args.dst_port)))
    if (float(args.udp_speed) > 0):
        time.sleep((file_duration - 1) / float(args.udp_speed))
    else:
        time.sleep(3600 * 4)
    raw_object.stop()
