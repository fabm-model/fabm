#!/usr/bin/env python

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
    parser = argparse.ArgumentParser(description='Convert Mauna Loa CO2 data to GOTM format')
    parser.add_argument('obs', type=str, help='Data file downloaded from: ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt')
    args = parser.parse_args()
    read_file(args.obs)

if __name__ == "__main__":
    main()

