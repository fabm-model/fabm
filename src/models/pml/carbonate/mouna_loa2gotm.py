#!/usr/bin/env python

def get_file(fn):
    from ftplib import FTP

    ftp = FTP('aftp.cmdl.noaa.gov')
    ftp.login()
    ftp.cwd('products/trends/co2')
    ftp.retrbinary('RETR co2_mm_mlo.txt', open(fn , 'wb').write)
    return

def read_file(fn):
    with open(fn,'rU') as f:
        while 1:
            l = f.readline()
            if l.startswith("#"):
                continue
            if l.strip()=='': break
            x = l.strip().split()
            y = int(x[0])
            m = int(x[1])
            v = float(x[4])
            print("%04d-%02d-15 00:00:00 %7.2f" % (y,m,v))

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert Mauna Loa CO2 data to GOTM format. Data are downloaded automatically via ftp.")
    args = parser.parse_args()
    fn = 'co2_mm_mlo.txt'
    get_file(fn)
    read_file(fn)

if __name__ == "__main__":
    main()

