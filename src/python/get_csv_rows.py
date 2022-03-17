#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import csv

debug = False

def get_csv_rows(csv_file, minimum_nr_columns, enc='utf-8'):
    delim     = "|"
    inputfile = open(csv_file, "r", encoding=enc)
    csv.register_dialect('pipes', delimiter=delim)

    rows = []
    try:
        reader = csv.reader(inputfile, dialect='pipes')

        for row in reader:
            if (len(row) >= minimum_nr_columns) and (row[0].find("#") < 0):
                rows.append(row)
    finally:
        inputfile.close()
    if debug:
        print("Read %d rows from %s" % ( len(rows), csv_file ))
    return rows

