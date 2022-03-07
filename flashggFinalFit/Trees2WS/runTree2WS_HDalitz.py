import os, sys

os.system("python tree2ws_HDalitz.py -ic config_HDalitz.py --year 2017 --mass 120")
os.system("python tree2ws_HDalitz.py -ic config_HDalitz.py --year 2017 --mass 125")
os.system("python tree2ws_HDalitz.py -ic config_HDalitz.py --year 2017 --mass 130")
os.system("python tree2ws_HDalitz_data.py -ic config_HDalitz_data.py")