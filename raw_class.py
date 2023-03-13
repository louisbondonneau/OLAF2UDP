
from pathlib import Path
import numpy as np
import re
import os
import sys
import glob

from astropy.time import Time, TimeDelta

from waveolaf_structure import WaveOlaf

import time
import threading

import socket

iframe = 0
ts_0 = 0
bsn_0 = 0
timestamp_0 = 0


class Raw:
    def __init__(self, obs_files, start=None, end=None, duration=None, block_start=0, block_end=1, verbose=False):
        self.block_start, self.block_end = block_start, block_end
        if(self.block_end < self.block_start):
            self.block_end, self.block_start = self.block_start, self.block_end
        self.names = obs_files
        self.fnames_raw = sorted(list(obs_files))
        self.data = []
        self.header = []
        self.verbose = verbose
        print(self.fnames_raw)

        fn_raw = Path(self.fnames_raw[0])
        if (WaveOlaf.try_waveolaf(fn_raw)):
            file_format = WaveOlaf(fn_raw, verbose=self.verbose)
        else:
            raise NameError('444: No suitable format found')

        self.file_format = file_format

        if(start is not None):
            if (self.__try_iso(start)):
                start = self.__try_iso(start)
            elif (self.__try_isot(start)):
                start = self.__try_isot(start)
            else:
                start = float(start)
        else:
            start = 0.0

        if(end is not None):
            if (self.__try_iso(end)):
                end = self.__try_iso(end)
            elif (self.__try_isot(end)):
                end = self.__try_isot(end)
            else:
                end = float(end)
        elif(duration is not None):
            end = start + duration
        else:
            end = 1.0

        if(start > end):
            end, start = start, end

        for file_num, (fn_raw) in enumerate(self.fnames_raw):
            # print(file_num)
            print(fn_raw)
            fn_raw = Path(fn_raw)
            if verbose:
                print("Opening %s" % (fn_raw.name))

            # dt_header, header, dt_block = self.infer_data_structure(fn_raw)
            data = file_format.data

            self.data.append(data)  # keep it if multiple file in future
            self.header.append(file_format.header)

            # header = self.header[file_num]

            if (file_num == 0):  # if first file
                self.nchan_file = file_format.nchan_file
                self.dm = file_format.dm
                self.nof_polcpx = file_format.nof_polcpx
                self.chan_bw = file_format.chan_bw
                self.imjd = file_format.imjd
                self.smjd = file_format.smjd
                self.offmjd = file_format.offmjd
                self.date_time = file_format.date_time

    def start(self, speed=1.0, adresse=("127.0.0.1", 20001)):
        self.adresse = adresse
        if (speed > 0):
            self.periodic_timer = PeriodicTimer(self.UDP_send_block, 5.12e-6 * float(self.file_format.nof_blocks_per_frame) / speed)
        else:
            self.periodic_timer = NoTimer(self.UDP_send_block, 5.12e-6)

    def stop(self):
        self.periodic_timer.stop()

    def UDP_send_block(self):
        global iframe
        if (iframe == 0):
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            # self.sock.connect(self.adresse)
            self.sock.setblocking(True)         # Now setting to non-blocking mode False in real time
        # if (iframe % 24415 == 0):
        nblock_per_sec = int(1 / self.file_format.tbin / self.file_format.nof_blocks_per_frame)
        if (iframe % nblock_per_sec == 0):
            print("iframe: %d -> %d sec" % (iframe, int(np.round(iframe * self.file_format.tbin * self.file_format.nof_blocks_per_frame))))
        # data = self.file_format.get_transposed_block(iframe)
        data = self.file_format.get_block(iframe)

        self.check_ts_bsn(data)
        # assert self.sock.send(data)
        assert self.sock.sendto(data, self.adresse)
        # ak_msg = self.sock.recvfrom(32)
        # print(ak_msg)
        iframe += 1

    def check_ts_bsn(self, data):
        global iframe, ts_0, bsn_0, timestamp_0

        header_tmp = np.frombuffer(data,
                                   count=1,
                                   dtype=self.file_format.dt_block,
                                   )[0]

        if (iframe == 0):
            ts_0 = int(header_tmp['Timestamp'])
            bsn_0 = int(header_tmp['BSN'])
            timestamp_0 = np.longdouble(ts_0) + (np.longdouble(bsn_0) * 5.12e-6)
        else:
            ts_1 = int(header_tmp['Timestamp'])
            bsn_1 = int(header_tmp['BSN'])
            nof_blocks_per_frame = int(header_tmp['nofblock'])
            timestamp_1 = np.longdouble(ts_1) + (np.longdouble(bsn_1) * 5.12e-6)
            if ((timestamp_1 - timestamp_0) > ((nof_blocks_per_frame + 1) * 5.12e-6)):
                print('Timestamp: WARNING deta dt = %.8f' % (timestamp_1 - (nof_blocks_per_frame * 5.12e-6) - timestamp_0))
                print('Timestamp: WARNING  ts_0: %d  ts_1: %d bsn_0: %d bsn_1: %d' % (ts_0, ts_1, bsn_0, bsn_1))
                missing_frame = int(np.round((timestamp_1 - timestamp_0) / (nof_blocks_per_frame * 5.12e-6))) - 1
                print('Missing frame: missing %d frame (%f sec) packets at %f sec' %
                      (missing_frame, missing_frame * nof_blocks_per_frame * 5.12e-6, iframe * nof_blocks_per_frame * 5.12e-6))
            timestamp_0 = timestamp_1
            ts_0 = ts_1
            bsn_0 = bsn_1

    def get_starttime(self):
        return self.file_format.date_time

    def get_stoptime(self):
        duration = self.file_format.nof_frame * self.file_format.nof_blocks_per_frame * self.file_format.tbin
        stop_time = self.get_starttime() + TimeDelta(float(duration), format='sec')
        return stop_time


class PeriodicTimer(threading.Thread):
    def __init__(self, task_function, period):
        super().__init__()
        self.task_function = task_function
        self.period = period
        self.i = 0
        self.t0 = time.time()
        self.RUN = True
        self.start()

    def __sleep(self):
        self.i += 1
        delta = self.t0 + self.period * self.i - time.time()
        if delta > 0:
            time.sleep(delta)

    def run(self):
        # time.sleep(4)
        while self.RUN:
            self.task_function()
            self.__sleep()

    def stop(self):
        self.RUN = False


class NoTimer(threading.Thread):
    def __init__(self, task_function, period):
        super().__init__()
        self.task_function = task_function
        self.period = period
        self.i = 0
        self.RUN = True
        self.start()

    def run(self):
        time.sleep(4)
        while self.RUN:
            self.task_function()

    def stop(self):
        self.RUN = False
