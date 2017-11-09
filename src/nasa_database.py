import sqlite3
import re

# create database
db = sqlite3.connect('NASAPolyCoeff.sqlite')
cursor = db.cursor()
cursor.execute("DROP TABLE IF EXISTS PolyCoeff")
cursor.execute("PRAGMA foreign_keys=1")
cursor.execute('''CREATE TABLE PolyCoeff (
               id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
               specie TEXT NOT NULL, 
               coeff TEXT NOT NULL, 
               temp TEXT NOT NULL)''')
db.commit()

# read from the text file
with open('thermo30.txt') as f:
    thermo = [item.strip() for item in f.readlines()]

thermo = thermo[5:-5]

# create a list from the data
thermo_list = []

for i, item in enumerate(thermo):
    if i%4 == 0:
        items = item.split()
        specie = items[0]
    elif i%4 == 1:
        items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
        high_temp = items
    elif i%4 == 2:
        items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
        high_temp += items[:2]
        low_temp = items[2:]
    else:
        items = re.findall('[-]?[0-9]+\.[0-9]+[E][+-][0-9]+', item)
        low_temp += items
        thermo_list.append([specie, ', '.join(low_temp), 'low'])
        thermo_list.append([specie, ', '.join(high_temp), 'high'])

# insert into the database
for i, item in enumerate(thermo_list):
    cursor.execute('''INSERT INTO PolyCoeff (specie, coeff, temp) 
                               VALUES (?, ?, ?)''', item)

db.commit()
db.close()