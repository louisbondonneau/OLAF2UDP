
import numpy as np
import struct
from datetime import datetime
from astropy.time import Time
from multiprocessing import Pool


def transpose_block(block):
    return block.transpose((1, 0, 2))


class WaveOlaf:
    def __init__(self, fn_raw, verbose=False):
        """Parse raw file
        """
        self.verbose = verbose
        self.header, self.dt_header = self.__header(fn_raw)

        self.npol = 4
        overlap = 0
        self.nchan = int(self.header['nof_beamlets'])
        self.nbit = int(8)
        self.nschan_file = int(self.header['nofblock']) * self.npol
        self.blocksize = self.nchan * self.nschan_file

        if (self.nbit == 8):
            self.bit_mode = 'uint8'
        elif (self.nbit == 16):
            self.bit_mode = 'uint16'
        else:
            print("Unknown nbit:")

        self.dt_block = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32'),
                                  ('data', self.bit_mode, (int(self.blocksize * self.nbit / 8))),
                                  ])

        self.dt_block_bytes = np.dtype([
            ('block', 'bytes', (self.dt_block.itemsize)),
        ])

        if self.verbose:
            print("  blocksize_header  is: %s -> %.8f sec" % (self.blocksize, float((self.blocksize / float(self.nchan) / self.npol) - overlap) * 5.12e-6))
            print("  blocksize_file    is: %s" % self.dt_block['data'].itemsize)
            print("  nbit       is: %s" % self.nbit)
            print("  head_size  is: %s" % self.dt_header.itemsize)

        self.data = self.__open_raw(fn_raw)
        self.nchan_file = int(self.header['nof_beamlets'])
        self.dm = 0
        self.nof_frame = int(self.data['block'].shape[0])
        self.nof_blocks_per_frame = int(self.header['nofblock'])
        self.nof_polcpx = int(4)
        self.overlap = 0
        self.tbin = float(5.12e-6)
        self.chan_bw = float(200. / 1024)

        date_time = datetime.fromtimestamp(self.header['Timestamp']).isoformat()
        self.date_time = Time(date_time, format='isot', scale='utc')

        self.imjd = np.floor(self.date_time.mjd)
        self.smjd = np.round((self.date_time.mjd % 1) * 3600 * 24)
        if(self.header['Timestamp'] % 1 != 0):
            self.offmjd = self.header['BSN'] * self.tbin + self.tbin / 2
        else:
            self.offmjd = self.header['BSN'] * self.tbin

    def get_patidx(self, file_num, block_num):
        ts_0 = self.data['Timestamp'][0]
        if(ts_0 % 1 != 0):
            ts_0 += self.tbin / 2
        ts_1 = self.data['Timestamp'][block_num]
        if(ts_1 % 1 != 0):
            ts_1 += self.tbin / 2

        ts_from_start = ts_1 - ts_0
        return int(ts_from_start / self.tbin + self.data['BSN'][block_num])

    def get_transposed_block2(self, block_num):
        unpack_block = np.frombuffer(self.get_block(block_num), dtype=self.dt_block)
        data = unpack_block["data"][0]
        data = data.reshape(int(self.nchan_file), int(self.nof_blocks_per_frame), self.nof_polcpx)
        pool = Pool(8)
        block_size = int(self.nof_blocks_per_frame * self.nchan_file * self.nof_polcpx / pool._processes)
        result = pool.map(transpose_block, [data[i:i + block_size] for i in range(0, len(data), block_size)])
        data = np.concatenate(result)
        unpack_block["data"][0] = data.ravel()
        block_bytes = struct.pack(f"<BHBHBBii{int(self.blocksize * self.nbit / 8)}B", *unpack_block)
        return block_bytes

    def get_transposed_block(self, block_num):
        # self.data['nof_beamlets'][block_num] = nchan
        unpack_block = np.frombuffer(self.get_block(block_num), dtype=self.dt_block)
        self.nbit
        self.blocksize
        struct_string = "<BHBHBBii%dB" % int(self.blocksize * self.nbit / 8)
        data = unpack_block["data"][0]  # size = self.nof_blocks_per_frame * self.nchan * self.npol
        data.shape = (int(self.nchan_file), int(self.nof_blocks_per_frame), self.nof_polcpx)
        data = np.ascontiguousarray(data.transpose((1, 0, 2)))
        data.shape = (int(self.nof_blocks_per_frame) * int(self.nchan_file) * self.nof_polcpx)

        block_bytes = struct.pack(struct_string, unpack_block["version_id"][0],
                                  unpack_block["source_info"][0],
                                  unpack_block["configuration_id"][0],
                                  unpack_block["station_id"][0],
                                  unpack_block["nof_beamlets"][0],
                                  unpack_block["nofblock"][0],
                                  unpack_block["Timestamp"][0],
                                  unpack_block["BSN"][0],
                                  *tuple(data))

        # out = self.data[block_num]
        # out = out.transpose((1, 0, 2))
        # return np.ascontiguousarray(out)  # shape is 1D (nchan, nschan, npol)
        return block_bytes  # shape is 1D (nchan, nschan, npol)

    def get_block(self, block_num):
        out = self.data[block_num]
        # out.shape = (int(self.nschan_file / self.nof_polcpx), self.nchan_file, self.nof_polcpx)
        # out = out.transpose((1, 0, 2))
        return out  # np.ascontiguousarray(out) #shape is 1D (nchan, nschan, npol)

    def __open_raw(self, fn_raw):
        with fn_raw.open('rb') as fd_raw:
            fd_raw.seek(0, 2)
            file_size = fd_raw.tell()
            nblock = int(np.floor(file_size / self.dt_block.itemsize))
            print("  nblock     is: %s" % nblock)
            if((fd_raw.tell() % float(self.dt_block.itemsize)) != 0):
                print('WARNING: last block is corrupted %f' % (fd_raw.tell() / float(self.dt_block.itemsize)))
        data = np.memmap(fn_raw.as_posix(),
                         dtype=self.dt_block_bytes,
                         # dtype=self.dt_block,
                         mode='r',
                         shape=(nblock,),
                         offset=0,
                         )
        return data

    def __header(self, fn_raw):
        dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32'),
                                  ])
        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]

        return header_tmp, dt_header_tmp

    def try_waveolaf(fn_raw):

        dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32')
                                  ])
        try:
            with fn_raw.open('rb') as fd_raw:
                header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                           count=1,
                                           dtype=dt_header_tmp,
                                           )[0]
            print("  nofblock per frame: %s" % header_tmp['nofblock'])
            print("  nof_beamlets: %s" % header_tmp['nof_beamlets'])
            blocksize = 4 * header_tmp['nofblock'] * header_tmp['nof_beamlets']
            dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                      ('source_info', 'uint16'),
                                      ('configuration_id', 'uint8'),
                                      ('station_id', 'uint16'),
                                      ('nof_beamlets', 'uint8'),
                                      ('nofblock', 'uint8'),
                                      ('Timestamp', 'int32'),
                                      ('BSN', 'int32'),
                                      ('data', 'uint8', (int(blocksize))),
                                      ])
            with fn_raw.open('rb') as fd_raw:
                fd_raw.seek(0, 2)
                file_size = fd_raw.tell()
                nblock = int(np.floor(file_size / dt_header_tmp.itemsize))
                print("  nblock     is: %s" % nblock)
                if((fd_raw.tell() % float(dt_header_tmp.itemsize)) != 0):
                    print('WARNING: last block is corrupted %f' % (fd_raw.tell() / float(dt_header_tmp.itemsize)))
            np.memmap(fn_raw.as_posix(),
                      dtype=dt_header_tmp,
                      mode='r',
                      offset=0,
                      shape=(nblock,)
                      )
        except ValueError:
            print("waveolaf is not a valid format.  Try again...")
            return False
        return True
