#!/usr/bin/env python3
"""
This script created the nasa polynomial coefficient database.

"""

import sqlite3
import re


if __name__ == "__main__":
    # create database
    db = sqlite3.connect('NASA_polynomial_coefficients.sqlite')
    cursor = db.cursor()
    cursor.execute("DROP TABLE IF EXISTS coefficients")
    cursor.execute('''CREATE TABLE coefficients (
                   id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                   species TEXT NOT NULL,
                   low_temp DECIMAL(10,2) NOT NULL,
                   mid_temp DECIMAL(10,2) NOT NULL,
                   high_temp DECIMAL(10,2) NOT NULL,
                   low_coeffs TEXT NOT NULL,
                   high_coeffs TEXT NOT NULL
                   )''')
    db.commit()

    # read from the text file
    with open('../src/thermo30.txt') as f:
        thermo = [item.strip() for item in f.readlines()]

    thermo = thermo[5:-5]

    # create a list from the data
    thermo_list = []

    for i, item in enumerate(thermo):
        if i % 4 == 0:
            low_t = float(item[45:55])
            high_t = float(item[55:65])
            mid_t = float(item[65:73])
            items = item.split()
            species = items[0]
        elif i % 4 == 1:
            items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
            high_temp = items
        elif i % 4 == 2:
            items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
            high_temp += items[:2]
            low_temp = items[2:]
        else:
            items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
            low_temp += items
            thermo_list.append(
                [species, low_t, mid_t, high_t,
                 ' '.join(low_temp), ' '.join(high_temp)])

    # insert into the database
    cursor.executemany('''
        INSERT INTO coefficients
        (species, low_temp, mid_temp, high_temp, low_coeffs, high_coeffs)
        VALUES (?, ?, ?, ?, ?, ?)''', thermo_list)

    db.commit()
    db.close()
